"""Example of using pyguppy_client_lib with read until"""
import argparse
import logging
import sys
import time

import numpy as np

from read_until import AccumulatingCache, ReadUntilClient

try:
    from pyguppy_client_lib.pyclient import PyGuppyClient
    from pyguppy_client_lib.helper_functions import package_read
except ImportError:
    print(
        "Failed to import pyguppy_client_lib, do you need to `pip install ont-pyguppy-client-lib`",
        file=sys.stderr,
    )


def basecall(
    guppy_client: PyGuppyClient, reads: list, dtype: "np.dtype", daq_values: dict,
):
    """Generator that sends and receives data from guppy

    :param guppy_client: pyguppy_client_lib.pyclient.PyGuppyClient
    :param reads: List of reads from read_until
    :type reads: Iterable
    :param dtype:
    :param daq_values:

    :returns:
        - read_info (:py:class:`tuple`) - channel (int), read number (int)
        - read_data (:py:class:`dict`) - Data returned from Guppy
    :rtype: Iterator[tuple[tuple, dict]]
    """
    hold = {}
    missing_count = 0

    with guppy_client:
        for channel, read in reads:
            hold[read.id] = (channel, read.number)
            t0 = time.time()
            success = guppy_client.pass_read(
                package_read(
                    read_id=read.id,
                    raw_data=np.frombuffer(read.raw_data, dtype),
                    daq_offset=daq_values[channel].offset,
                    daq_scaling=daq_values[channel].scaling,
                )
            )
            if not success:
                logging.warning("Skipped a read: {}".format(read.id))
                hold.pop(read.id)
                continue
            else:
                missing_count += 1

            sleep_time = guppy_client.throttle - t0
            if sleep_time > 0:
                time.sleep(sleep_time)

        while missing_count:
            results = guppy_client.get_completed_reads()
            missing_count -= len(results)

            if not results:
                time.sleep(guppy_client.throttle)
                continue

            yield from iter(
                [(hold[read["metadata"]["read_id"]], read) for read in results]
            )


def get_parser():
    """Build argument parser for example"""
    parser = argparse.ArgumentParser(prog="demo ({})".format(__file__),)
    parser.add_argument("--alignment-index-file", type=str, required=True)
    parser.add_argument("--bed-file", type=str, default="")
    parser.add_argument("--barcodes", nargs="*", default=[], help="Barcode kits in use")

    parser.add_argument(
        "--host", default="127.0.0.1", help="MinKNOW server host address"
    )
    parser.add_argument(
        "--port", type=int, default=8000, help="MinKNOW gRPC server port"
    )

    parser.add_argument(
        "--guppy_host", default="127.0.0.1", help="Guppy server host address",
    )
    parser.add_argument(
        "--guppy_port", type=int, default=5555, help="Guppy server port",
    )
    parser.add_argument(
        "--guppy_config", default="dna_r9.4.1_450bps_fast", help="Guppy server config",
    )

    parser.add_argument(
        "--run_time", type=int, default=30, help="Period to run the analysis"
    )
    parser.add_argument(
        "--unblock_duration",
        type=float,
        default=0.1,
        help="Time (in seconds) to apply unblock voltage",
    )
    parser.add_argument(
        "--batch_size",
        default=None,
        type=int,
        help="Number of reads to get from ReadCache each iteration. If not set uses number of channels on device",
    )
    parser.add_argument(
        "--throttle",
        default=0.1,
        type=float,
        help="Time to wait between requesting successive read batches from the ReadCache",
    )
    return parser


def analysis(
    client: ReadUntilClient,
    caller: PyGuppyClient,
    batch_size: int,
    duration: int,
    throttle: float = 0.1,
    unblock_duration: float = 0.1,
):
    """Example analysis function

    This is an example analysis function that collects data (read chunks)
    from the `ReadUntilClient`, passes them to the `PyGuppyClient` for
    calling, aligning, or barcoding and iterates the results.

    :param client: an instance of a `ReadUntilClient` object.
    :param caller: PyGuppyClient
    :param batch_size: number of reads to pull from `client` at a time.
    :param duration: time to run for, seconds
    :param throttle: minimum interval between requests to `client`.
    :param unblock_duration: time in seconds to apply unblock voltage.
    """
    run_duration = time.time() + duration
    logger = logging.getLogger("Analysis")
    # Get whether a bed file is in use from the pyguppy caller

    while client.is_running and time.time() < run_duration:
        time_begin = time.time()
        # Get most recent read chunks from read until client
        read_batch = client.get_read_chunks(batch_size=batch_size, last=True)

        # Send read_batch for base calling
        called_batch = basecall(
            guppy_client=caller,
            reads=read_batch,
            dtype=client.signal_dtype,
            daq_values=client.calibration_values,
        )

        n = 0
        for (channel, read_number), read in called_batch:
            n += 1
            alignment_genome = read["metadata"].get("alignment_genome")
            bed_hits = read["metadata"].get("bed_hits")
            barcode = read["metadata"].get("barcode_arrangement")
            barcode_kit = read["metadata"].get("barcode_kit")
            logger.debug(
                "{},{},{},{},{}".format(
                    read["metadata"]["read_id"],
                    alignment_genome,
                    bed_hits,
                    barcode,
                    barcode_kit,
                )
            )
            client.stop_receiving_read(channel, read_number)
            client.unblock_read(channel, read_number, unblock_duration)

        batch_time = time.time() - time_begin
        if n:
            logger.info("Processed {} in {:.5f}s".format(n, batch_time))

        # Limit the rate at which we make requests
        if batch_time < throttle:
            time.sleep(throttle - batch_time)


def main(argv=None):
    """simple example main cli entrypoint"""
    args = get_parser().parse_args(argv)

    logging.basicConfig(
        format="[%(asctime)s - %(name)s] %(message)s", level=logging.INFO,
    )

    read_until_client = ReadUntilClient(
        mk_host=args.host,
        mk_port=args.port,
        cache_type=AccumulatingCache,
        one_chunk=False,
        filter_strands=True,
        # Request uncalibrated, int16, signal
        calibrated_signal=False,
    )

    # Handle arg cases:
    if args.batch_size is None:
        args.batch_size = read_until_client.channel_count

    caller = PyGuppyClient(
        address="{}:{}".format(args.guppy_host, args.guppy_port),
        config=args.guppy_config,
        alignment_index_file=args.alignment_index_file,
        bed_file=args.bed_file if args.bed_file else "",
        barcode_kits=args.barcodes,
        server_file_load_timeout=180,  # 180 == 3 minutes, should be enough?
    )

    caller.connect()
    read_until_client.run()

    try:
        analysis(
            client=read_until_client,
            caller=caller,
            duration=args.run_time,
            batch_size=args.batch_size,
            unblock_duration=args.unblock_duration,
        )
    except KeyboardInterrupt:
        pass
    finally:
        read_until_client.reset()
        caller.disconnect()


if __name__ == "__main__":
    main()
