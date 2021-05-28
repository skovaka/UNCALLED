"""
MinKNOW RPC Access
==================

Provides access to MinKNOW via RPC.

This RPC system is gRPC-based. You might want to look at the `gRPC documentation
<http://www.grpc.io/grpc/python/>`_ for more information, but most of the detail is hidden by the
code in this module.

For external systems accessing MinKNOW, start with `minknow_api.manager.Manager` - this will allow
you to access high-level information, and enumerate and access each flow cell position.

The central class for accessing a flow cell position is `Connection`.

The functionality available on a flow cell position is divided into services, covering related units
like controlling the device or managing data acquisition. Each service is available as a property on
a `Connection` object. For each service, the related Protobuf messages are available from
``minknow_api.<service>_service``, or as ``connection.<service>._pb`` (if ``connection`` is a
``Connection`` object).

.. _rpc-services:

Services
--------

The available services are:

acquisition
    Control data acquisition. See `acquisition_service.AcquisitionService` for a description of the
    available methods.
analysis_configuration
    Configure data acquisition. See `analysis_configuration_service.AnalysisConfigurationService`
    for a description of the available methods.
data
    Stream acquisition data. Note that this is for data directly produced during acquisition, rather
    than statistics about acquired data. See `data_service.DataService` for a description of the
    available methods.
device
    Get information about and control the attached device. This useful presents information and
    settings in a device-independent way, so it can be used on PromethIONs as easily as on MinIONs.
    See `device_service.DeviceService` for a description of the available methods.
keystore
    A service for storing and retreiving arbitrary data on the instance. This can be used to
    communicate with other users of the API. See `keystore_service.DeviceService` for a description
    of the available methods.
instance
    Get information about the instance of MinKNOW you are connected to (eg: software version). See
    `instance_service.InstanceService` for a description of the available methods.
log
    Get or produce general informational messages. See `log_service.LogService` for a description of
    the available methods.
minion_device
    MinION-specific device interface. This exposes low-level settings for MinIONs and similar
    devices (eg: GridIONs). See `minion_device_service.MinionDeviceService` for a
    description of the available methods.
protocol
    Control protocol scripts. See `protocol_service.ProtocolService` for a description of the
    available methods.
promethion_device
    PromethION-specific device interface. This exposes low-level settings for PromethIONs. See
    `minion_device_service.MinionDeviceService` for a description of the available methods.
statistics
    Get statistics about an acquisition period. Statistics can be streamed live during acquisition,
    or retreived afterwards. See `statistics_service.StatisticsService` for a description of the
    available methods.


Helpers
-------

The `minknow_api.data`, `minknow_api.device` and `minknow_api.manager` modules contain helpers for
working with the services of the same name. See the documentation for those modules for more
information.

"""

import grpc
import logging
import os

#
# Services
#
_services = {
    "acquisition": ["AcquisitionService"],
    "analysis_configuration": ["AnalysisConfigurationService"],
    "data": ["DataService"],
    "device": ["DeviceService"],
    "keystore": ["KeyStoreService"],
    "instance": ["InstanceService"],
    "log": ["LogService"],
    "minion_device": ["MinionDeviceService"],
    "protocol": ["ProtocolService"],
    "production": ["ProductionService"],
    "promethion_device": ["PromethionDeviceService"],
    "statistics": ["StatisticsService"],
}
_optional_services = ["production"]


#
# Module meta-information
#

__all__ = [svc + "_service" for svc in _services] + [
    "Connection",
    "data",
    "device",
    "grpc_credentials",
    "manager",
]

try:
    from ._version import __version__
except ImportError:
    __version__ = None


#
# Submodule imports
#

# Convenience imports for each service
import importlib

for svc in _services:
    try:

        # effectively does `import .{svc}_service as {svc}_service`
        importlib.import_module(".{}_service".format(svc), __name__)
    except ImportError:
        if svc not in _optional_services:
            raise

from . import data
from . import manager

logger = logging.getLogger(__name__)

#
# Connection helpers
#


class MissingMinknowSSlCertError(Exception):
    pass


def ca_path_from_minknow():
    try:
        from minknow.paths import minknow_base_dir
    except ImportError:
        return None

    return os.path.join(minknow_base_dir(), "conf", "rpc-certs", "ca.crt")


def get_ca_path():
    try:
        cert_path = os.environ["MINKNOW_TRUSTED_CA"]
    except KeyError:
        cert_path = ca_path_from_minknow()

    return cert_path


def grpc_credentials():
    """Get a grpc.ChannelCredentials object for connecting to secure versions of MinKNOW"s gRPC
    services.

    Use like:

    >>> import grpc
    >>> channel = grpc.secure_channel("localhost:9502", grpc_credentials())

    If run from the Python embedded in MinKNOW, this will find the correct certificate
    automatically. Otherwise, you may need to set the ``MINKNOW_TRUSTED_CA`` to point to the CA
    certificate used by MinKNOW (which can be found at ``conf/rpc-certs/ca.crt`` in the MinKNOW
    installation).
    """
    try:
        return grpc_credentials.cached_credentials
    except AttributeError:
        cert_path = get_ca_path()

        if not cert_path:
            raise MissingMinknowSSlCertError(
                "Couldn't find a valid path to MinKNOW's CA SSL certificate to initiate a secure connection"
            )

        # cert_path = os.path.join(minknow.paths.minknow_base_dir(), "conf", "rpc-certs", "ca.crt")
        with open(cert_path, "rb") as file_obj:
            grpc_credentials.cached_credentials = grpc.ssl_channel_credentials(
                file_obj.read()
            )
        return grpc_credentials.cached_credentials


#
# Connection class
#


class Connection(object):
    """A connection to a MinKNOW sequencing position via RPC.

    Note that this only provides access to the new RPC system. The old one is
    available via the minknow.engine_client module.

    Each service is available as a property of the same name on the Connection
    object. See :ref:`rpc-services` for a list.

    Given a connection object ``connection``, for each service,
    ``connection.<service>`` is a "service object". This exposes the RPC methods
    for that service in a more convenient form than gRPC's own Python bindings
    do.

    For example, when calling ``start_protocol`` on the ``protocol`` service,
    instead of doing

    >>> protocol_service.start_protocol(
    >>>     protocol_service._pb.StartProtocolMessage(path="my_script"))

    you can do

    >>> connection.protocol.start_protocol(path="my_script")

    Note that you must use keyword arguments - no positional parameters are
    available.

    This "unwrapping" of request messages only happens at one level, however. If
    you want to change the target temperature settings on a MinION, you need to do
    something like

    >>> temp_settings = connection.minion_device._pb.TemperatureRange(min=37.0, max=37.0)
    >>> connection.minion_device.change_settings(
    >>>     temperature_target=temp_settings)
    """

    def __init__(self, port=None, host="127.0.0.1", use_tls=False):
        """Connect to MinKNOW.

        The port for a given instance of MinKNOW is provided by the manager
        service.

        If no port is provided, it will attempt to get the port from the
        MINKNOW_RPC_PORT environment variable (set for protocol scripts, for
        example). If this environment variable does not exist (or is not a
        number), it will raise an error.

        :param port: the port to connect on (defaults to ``MINKNOW_RPC_PORT`` environment variable)
        :param host: the host MinKNOW is running on (defaults to localhost)
        """
        import grpc, os, time

        self.host = host
        if port is None:
            port = int(os.environ["MINKNOW_RPC_PORT"])
        self.port = port

        channel_opts = [
            ("grpc.max_send_message_length", 16 * 1024 * 1024),
            ("grpc.max_receive_message_length", 16 * 1024 * 1024),
            ("grpc.http2.min_time_between_pings_ms", 1000),
            (
                "grpc.ssl_target_name_override",
                "localhost",
            ),  # that's what our cert's CN is
        ]

        error = None
        retry_count = 5
        for i in range(retry_count):
            if use_tls:
                self.channel = grpc.secure_channel(
                    "{}:{}".format(host, port),
                    credentials=grpc_credentials(),
                    options=channel_opts,
                )
            else:
                self.channel = grpc.insecure_channel(
                    "{}:{}".format(host, port), options=channel_opts
                )

            # One entry for each service
            for name, svc_list in _services.items():
                for svc in svc_list:
                    try:
                        # effectively does `self.{name} = {name}_service.{svc}(self.channel)`
                        setattr(
                            self,
                            name,
                            getattr(globals()[name + "_service"], svc)(self.channel),
                        )
                    except KeyError:
                        if name not in _optional_services:
                            raise

            # Ensure channel is ready for communication
            try:
                self.instance.get_version_info()
                error = None
                break
            except grpc.RpcError as e:
                logger.info("Error received from rpc")
                if (
                    e.code() == grpc.StatusCode.INTERNAL
                    and e.details() == "GOAWAY received"
                ):
                    logger.warning(
                        "Failed to connect to minknow instance (retry %s/%s): %s",
                        i + 1,
                        retry_count,
                        e.details(),
                    )
                elif e.code() == grpc.StatusCode.UNAVAILABLE:
                    logger.warning(
                        "Failed to connect to minknow instance (retry %s/%s): %s",
                        i + 1,
                        retry_count,
                        e.details(),
                    )
                else:
                    raise
                error = e
                time.sleep(0.5)

        if error:
            raise error

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.channel.close()
