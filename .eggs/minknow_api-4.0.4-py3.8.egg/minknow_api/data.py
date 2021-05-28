"""
Helpers for accessing the data service
======================================

The primary function this module provides is `get_signal` for getting live signal data from MinKNOW.
The ``get_signal_bytes`` RPC has been designed to be performant to use from Python, but this also
makes it complex to use properly. `get_signal` provides an easy-to-use wrapper for the common case
where you just want to obtain a fixed amount of data, to be processed after you have it all.
"""

import collections
import numpy

from ._support import ArgumentError

__all__ = [
    "ChannelConfigChange",
    "ChannelSignalData",
    "NumpyDTypes",
    "SignalData",
    "api_types_to_numpy_types",
    "get_numpy_types",
    "get_signal",
]

ChannelConfigChange = collections.namedtuple(
    "ChannelConfigChange", ["offset", "config"]
)
ChannelConfigChange.__doc__ = """\
A channel configuration change for data in a ChannelSignalData.

Attributes:
(each of which has an ``offset`` field, which is an offset into the
        signal array, and a ``config`` field, which contains the updated configuration for the
        channel)
"""

ChannelSignalData = collections.namedtuple(
    "ChannelSignalData", ["name", "signal", "config_changes",]
)
ChannelSignalData.__doc__ = """\
The per-channel data in SignalData.

Attributes:
    name (int): The channel this is the data for (channel numbers start at 1).
    signal (numpy.ndarray): The signal data, as a 1-dimensional numpy array. If calibrated data was
        requested, this will be floating-point data in picoamps. Otherwise, this will be integer
        data (ADC values). Note that the exact type will depend on the sequencing device (and the
        host machine in the case of MinIONs - MinITs will use a slightly different format to PCs).
    config_changes (list of ChannelConfigChange): If channel configuration changes were requested, a
        list of those changes . This will contain at least one element, with offset 0, which is the configuration
        that applies to the first sample returned.
"""

NumpyDTypes = collections.namedtuple(
    "NumpyDTypes", ["bias_voltages", "calibrated_signal", "uncalibrated_signal",]
)
NumpyDTypes.__doc__ = """\
The types of data returned by the data.get_signal_bytes() RPC, as numpy dtypes.

See api_types_to_numpy_types() and get_numpy_types().

Note that there are no guarantees about the endianness of these data types.

Attributes:
    bias_voltages (numpy.dtype): The type of bias voltage data. This may be floating-point or
        integer type, depending on the type of sequencing device.
    calibrated_signal (numpy.dtype): The type of calibrated signal data. This will be a
        floating-point type.
    uncalibrated_signal (numpy.dtype): The type of uncalibrated (raw) signal data. This will be an
        integer type.
"""

SignalData = collections.namedtuple(
    "SignalData",
    ["samples_since_start", "seconds_since_start", "channels", "bias_voltages"],
)
SignalData.__doc__ = """\
The results of get_signal().

Attributes:
    samples_since_start (int): The number of samples collected before the first returned sample.
    seconds_since_start (int): As samples_since_start, but expressed in seconds.
    channels (list of ChannelSignalData): The per-channel data (signal and configuration changes).
    bias_voltages (numpy.ndarray): A numpy array of bias voltages in millivolts. This will
        be the same length as the ``signal`` array on each channel. Note that there should be no
        need to apply any further corrections to the value (eg: the 5x amplifier on a MinION is
        already accounted for). Be aware that the types stored in this array will be different for
        MinION-like devices (integers) and PromethION-like devices (floating-point).

        May be empty if no bias voltages were requested.
"""


def _numpy_type(desc):
    """Convert a description of a type provided by the data.get_data_types() RPC into a numpy dtype
    object."""

    if desc.type == desc.SIGNED_INTEGER:
        type_char = "i"
    elif desc.type == desc.UNSIGNED_INTEGER:
        type_char = "u"
    elif desc.type == desc.FLOATING_POINT:
        type_char = "f"
    else:
        raise RuntimeError("Unknown type format {}".format(desc))
    type_desc = "{}{}{}".format(">" if desc.big_endian else "<", type_char, desc.size)
    return numpy.dtype(type_desc)


def api_types_to_numpy_types(data_types):
    """Convert the result of data.get_data_types() into numpy types.

    Args:
        data_types (minknow_api.data_pb2.GetDataTypesResponse): The result of a call to the
            data.get_data_types() RPC.

    Result:
        NumpyDTypes: The numpy dtypes recorded in `data_types`.
    """
    return NumpyDTypes(
        _numpy_type(data_types.bias_voltages),
        _numpy_type(data_types.calibrated_signal),
        _numpy_type(data_types.uncalibrated_signal),
    )


def get_numpy_types(connection):
    """The data types provided by the data RPC service, in numpy format.

    Args:
        connection (minknow_api.Connection): Connection to a MinKNOW flow cell position.

    Result:
        NumpyDTypes: The numpy dtypes that will be returned over `connection`.
    """
    data_types = connection.data.get_data_types()
    return api_types_to_numpy_types(data_types)


def get_signal(connection, on_started=None, numpy_dtypes=None, **kwargs):
    """Get signal data from the flow cell.

    This can be used to sample the signal being produced by the flow cell. The signal can be
    returned as raw ADC values or as calibrated picoamp (pA) values; see ``set_calibration`` on
    the device service for the values used in this conversion.

    In addition to the signal, this can return the associated channel configuration and/or bias
    voltage information, to help analyse the data.

    If a device settings change RPC has completed before this method is called, the data returned
    is guaranteed to have been generated by the device after those settings were applied.
    However, note that no guarantee is made about how device settings changes that overlap with
    this request will affect the returned data.

    Exactly one of ``seconds`` or ``samples`` must be provided.

    Args:
        connection (minknow_api.Connection): Connection to a MinKNOW flow cell position.
        on_started (callable, optional): Called as soon as the first message is received.
        numpy_dtypes (NumpyDTypes, optional): The result of ``get_numpy_types(connection)``.
            Providing this avoids an extra RPC to determine this information - if making frequent
            calls to `get_signal`, it may be worth caching this information.
        seconds (float): Amount of data to collect in seconds.
        samples (int): Amount of data to collect in samples.
        first_channel (int): The first channel to collect data from. Channels start at 1.
        last_channel (int): The last channel to collect data from.
        include_channel_configs (bool): Whether to return changes in channel configurations.
        include_bias_voltages (bool): Whether to return bias voltage information.
        calibrated_data (bool): Return the data in picoamps rather than raw ADC values. Requires
            that a calibration has been set.
        return_when_listening (bool): Returns an empty first message as soon as the request is set
            up, possibly before there is data to return. Potentially useful in combination with
            `on_started`.

    Returns:
        SignalData: The returned data.
    """
    if numpy_dtypes is None:
        numpy_dtypes = get_numpy_types(connection)

    if "samples" not in kwargs and "seconds" not in kwargs:
        raise ArgumentError("Expected 'samples' or 'seconds' argument")

    if kwargs.get("calibrated_data", False):
        signal_dtype = numpy_dtypes.calibrated_signal
    else:
        signal_dtype = numpy_dtypes.uncalibrated_signal

    # don't raise here when these keys are missing - let the RPC call raise the correct error
    # instead
    first_channel = kwargs.get("first_channel", 0)
    channel_count = kwargs.get("last_channel", 0) + 1 - first_channel

    bias_voltages = []
    channel_configs = [[] for i in range(channel_count)]
    signal = [[] for i in range(channel_count)]

    start_samples = None

    for msg in connection.data.get_signal_bytes(**kwargs):

        if on_started:
            on_started()
            on_started = None

        if start_samples is None:
            offset = 0
            start_samples = msg.samples_since_start
            start_seconds = msg.seconds_since_start
        else:
            offset = msg.samples_since_start - start_samples
        for i, c in enumerate(msg.channels, start=msg.skipped_channels):
            signal[i].append(c.data)
            if len(c.config_changes):
                for change in c.config_changes:
                    channel_configs[i].append(
                        ChannelConfigChange(change.offset + offset, change.config)
                    )
        if len(msg.bias_voltages):
            bias_voltages.append(msg.bias_voltages)

    return SignalData(
        start_samples,
        start_seconds,
        [
            ChannelSignalData(
                channel, numpy.frombuffer(b"".join(ch_signal), signal_dtype), configs
            )
            for channel, (ch_signal, configs) in enumerate(
                zip(signal, channel_configs), start=first_channel
            )
        ],
        numpy.frombuffer(b"".join(bias_voltages), numpy_dtypes.bias_voltages),
    )
