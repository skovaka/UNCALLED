"""Tools for interacting with protocols in minknow"""

import collections
import logging
import typing

import grpc

from minknow_api import protocol_pb2

from .. import Connection

LOGGER = logging.getLogger(__name__)


def find_protocol(
    device_connection: Connection,
    product_code: str,
    kit: str,
    basecalling: bool = False,
    basecall_config: typing.Optional[str] = None,
    barcoding: bool = False,
    barcoding_kits: typing.Optional[typing.List[str]] = None,
    force_reload: bool = False,
    experiment_type: str = "sequencing",
) -> typing.Optional[str]:
    """Find a protocol identifier.

    This will fetch a list of protocols from the device-instance, then search through the protocols
    for one that supports the flow-cell type (product code) and all the specified options. It
    returns the first protocol it finds that matches.

    Args:
        device_connection (:obj:`Connection`):  As returned by minknow.manager.FlowCellPosition().connect().
        product_code (:obj:`str`):              The flow-cell type, as in flow_cell_info.product_code.
        kit (:obj:`str`):                       The kit to be sequenced. eg: "SQK-LSK108".
        basecalling (bool):                     True if base-calling is required
        basecall_config (:obj:`str):            The base-calling model that the protocol should support. If absent,
                                                the protocol default will be used (if specified)
        barcoding (bool):                       True if barcoding is required.
        barcoding_kits (:obj:`list`):           Barcoding kits that the protocol should support. If specified,
                                                barcoding is assumed to be True.        
        force_reload (bool):                    If true will force reloading the protocols from their descriptions,
                                                this will take a few seconds.
        experiment_type(:obj:`str`):            Type of experiment to be run.

    Returns:
        The first protocol to match or None.
    """

    try:
        response = device_connection.protocol.list_protocols(force_reload=force_reload)
    except grpc.RpcError as exception:
        raise Exception(
            "Could not get a list of protocols ({})".format(exception.details())
        )

    if not response.protocols:
        raise Exception("List of protocols is empty")

    for protocol in response.protocols:
        # we need the tags, if we don't have them move ont next protocol
        if not protocol.tag_extraction_result.success:
            LOGGER.debug("Ignoring protocol with tag extraction failure")
            continue

        # the tags provide a way to filter the experiments
        tags = protocol.tags

        # want a sequencing experiment...
        if experiment_type and tags["experiment type"].string_value != experiment_type:
            LOGGER.debug(
                "Ignoring experiment with incorrect type: %s vs %s",
                tags["experiment type"].string_value,
                experiment_type,
            )
            continue

        # ...for the correct flow-cell type...
        if tags["flow cell"].string_value != product_code:
            LOGGER.debug(
                "Protocol is for product %s, not %s",
                tags["flow cell"].string_value,
                product_code,
            )
            continue

        # ...with matching kit
        if tags["kit"].string_value != kit:
            LOGGER.debug(
                "Protocol supports kit %s, not %s", tags["kit"].string_value, kit
            )
            continue

        # if bar-coding is required, the protocol should support it and all
        # the bar-coding kits in use
        if tags["barcoding"].bool_value != barcoding:
            if barcoding:
                LOGGER.debug("Protocol does not support barcoding")
            else:
                LOGGER.debug("Protocol requires barcoding")
            continue

        supported_kits = tags["barcoding kits"].array_value
        # workaround for the set of barcoding kits being returned as a string rather
        # that array of strings
        if supported_kits and len(supported_kits[0]) == 1:
            supported_kits = (
                tags["barcoding kits"].array_value[1:-1].replace('"', "").split(",")
            )
        if barcoding_kits and not set(barcoding_kits).issubset(supported_kits):
            LOGGER.debug(
                "barcoding kits specified %s not amongst those supported %s",
                barcoding_kits,
                supported_kits,
            )
            continue

        # check base-calling is supported if required and with the requested model if
        # specified or default if not
        basecalling_required = basecalling or basecall_config is not None
        if basecalling_required:
            if not "default basecall model" in tags:
                continue
            default_basecall_model = tags.get("default basecall model").string_value
            if not basecall_config:
                basecall_config = default_basecall_model
            if basecall_config not in tags["available basecall models"].array_value:
                LOGGER.debug(
                    'basecalling model "%s" is not in the list of supported basecalling models %s',
                    basecall_config,
                    tags["available basecall models"].array_value,
                )
                continue

        # we have a match, (ignore the rest)
        return protocol
    return None


BarcodingArgs = collections.namedtuple(
    "BarcodingArgs",
    [
        "kits",
        "trim_barcodes",
        "barcodes_both_ends",
        "detect_mid_strand_barcodes",
        "min_score",
        "min_score_rear",
        "min_score_mid",
    ],
)
AlignmentArgs = collections.namedtuple("AlignmentArgs", ["reference_files", "bed_file"])
BasecallingArgs = collections.namedtuple(
    "BasecallingArgs", ["config", "barcoding", "alignment"]
)
OutputArgs = collections.namedtuple("OutputArgs", ["reads_per_file"])


def make_protocol_arguments(
    experiment_duration: float = 72,
    basecalling: BasecallingArgs = None,
    fastq_arguments: OutputArgs = None,
    fast5_arguments: OutputArgs = None,
    bam_arguments: OutputArgs = None,
    disable_active_channel_selection: bool = False,
    mux_scan_period: float = 1.5,
    args: typing.Optional[typing.List[str]] = None,
    is_flongle: bool = False,
) -> typing.List[str]:
    """Build arguments to be used when starting a protocol.
    
    This will assemble the arguments passed to this script into arguments to pass to the protocol.

    Args:
        experiment_duration(float):             Length of the experiment in hours.
        basecalling(:obj:`BasecallingArgs`):    Arguments to control basecalling.
        fastq_arguments(:obj:`OutputArgs`):     Control fastq file generation.
        fast5_arguments(:obj:`OutputArgs`):     Control fastq file generation.
        bam_arguments(:obj:`OutputArgs`):       Control bam file generation.
        disable_active_channel_selection(bool): Disable active channel selection
        mux_scan_period(float):                 Period of time between mux scans in hours.
        args(:obj:`list`):                      Extra arguments to pass to protocol.
        is_flongle(bool):                       Specify if the flow cell to be sequenced on is a flongle.

    Returns:
        A list of strings to be passed as arguments to start_protocol.
    """

    def on_off(value: bool):
        if value:
            return "on"
        else:
            return "off"

    protocol_args = []

    if basecalling:
        protocol_args.append("--base_calling=on")

        if basecalling.config:
            protocol_args.append("--guppy_filename=" + basecalling.config)

        if basecalling.barcoding:
            barcoding_args = []
            if basecalling.barcoding.kits:
                # list of barcoding kits converted to quoted, comma separated array elements
                # eg: barcoding_kits=['a','b','c']
                barcoding_args.append(
                    "barcoding_kits=['" + "','".join(basecalling.barcoding.kits) + "',]"
                )

            if basecalling.barcoding.trim_barcodes:
                # trim_barcodes=on/off
                barcoding_args.append(
                    "trim_barcodes=" + on_off(basecalling.barcoding.trim_barcodes)
                )

            if basecalling.barcoding.barcodes_both_ends:
                # require_barcodes_both_ends=on/off
                barcoding_args.append(
                    "require_barcodes_both_ends="
                    + on_off(basecalling.barcoding.barcodes_both_ends)
                )

            if basecalling.barcoding.detect_mid_strand_barcodes:
                # detect_mid_strand_barcodes=on/off
                barcoding_args.append(
                    "detect_mid_strand_barcodes="
                    + on_off(basecalling.barcoding.detect_mid_strand_barcodes)
                )

            if basecalling.barcoding.min_score:
                # min_score=66
                barcoding_args.append(
                    "min_score={}".format(basecalling.barcoding.min_score)
                )

            if basecalling.barcoding.min_score_rear:
                # min_score_rear=66
                barcoding_args.append(
                    "min_score_rear={}".format(basecalling.barcoding.min_score_rear)
                )

            if basecalling.barcoding.min_score_mid:
                # min_score_mid=66
                barcoding_args.append(
                    "min_score_mid={}".format(basecalling.barcoding.min_score_mid)
                )

            protocol_args.extend(["--barcoding"] + barcoding_args)

        if basecalling.alignment:
            alignment_args = []
            if basecalling.alignment.reference_files:
                alignment_args.append(
                    "reference_files=['"
                    + "','".join(basecalling.alignment.reference_files)
                    + "',]"
                )
            if basecalling.alignment.bed_file:
                alignment_args.append(
                    "bed_file='{}'".format(basecalling.alignment.bed_file)
                )
            protocol_args.extend(["--alignment"] + alignment_args)

    protocol_args.append("--experiment_time={}".format(experiment_duration))
    protocol_args.append("--fast5=" + on_off(fast5_arguments))
    if fast5_arguments:
        protocol_args.extend(
            ["--fast5_data", "trace_table", "fastq", "raw", "vbz_compress"]
        )
        protocol_args.append(
            "--fast5_reads_per_file={}".format(fast5_arguments.reads_per_file)
        )

    protocol_args.append("--fastq=" + on_off(fastq_arguments))
    if fastq_arguments:
        protocol_args.extend(["--fastq_data", "compress"])
        protocol_args.append(
            "--fastq_reads_per_file={}".format(fastq_arguments.reads_per_file)
        )

    protocol_args.append("--bam=" + on_off(bam_arguments))
    if bam_arguments:
        if bam_arguments.reads_per_file != 4000:
            raise Exception("Unable to change reads per file for BAM.")
        """protocol_args.append(
            "--bam_reads_per_file={}".format(bam_arguments.reads_per_file)
        )"""

    if not is_flongle:
        protocol_args.append(
            "--active_channel_selection=" + on_off(not disable_active_channel_selection)
        )
        if not disable_active_channel_selection:
            protocol_args.append("--mux-scan-period={}".format(mux_scan_period))

    protocol_args.extend(args)

    return protocol_args


def start_protocol(
    device_connection: Connection,
    identifier: str,
    sample_id: str,
    experiment_group: str,
    *args,
    **kwargs
) -> str:
    """Start a protocol on the passed {device_connection}.

    Args:
        device_connection(:obj:`Connection`):   The device connection to start a protocol on.
        identifier(str):                        Protocol identifier to be started.
        sample_id(str):                         Sample id of protocol to start.
        experiment_group(str):                  Experiment group of protocol to start.
        *args: Additional arguments forwarded to {make_protocol_arguments}
        **kwargs: Additional arguments forwarded to {make_protocol_arguments}

    Returns:
        The protocol_run_id of the started protocol.
    """

    flow_cell_info = device_connection.device.get_flow_cell_info()

    protocol_arguments = make_protocol_arguments(
        *args, is_flongle=flow_cell_info.has_adapter, **kwargs
    )
    LOGGER.debug("Built protocol arguments: %s", " ".join(protocol_arguments))

    user_info = protocol_pb2.ProtocolRunUserInfo()
    if sample_id:
        user_info.sample_id.value = sample_id
    if experiment_group:
        user_info.protocol_group_id.value = experiment_group

    result = device_connection.protocol.start_protocol(
        identifier=identifier, args=protocol_arguments, user_info=user_info
    )

    return result.run_id
