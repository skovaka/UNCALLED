"""
Helpers for accessing the manager
=================================

If you are connecting to MinKNOW externally (rather than from a protocol script), you will need to
go via the manager. This listens on a known port, and can enumerate the available flow cell
positions. There may only be one position (eg: on a Mk1C), or there may be many (PromethIONs can
have 24 or 48, for example).

This can be done with the `Manager` class in this module. The other classes this module provides are
usually constructing using methods on ``Manager``.

"""


import grpc
import logging
import minknow_api
import minknow_api.basecaller_service
import minknow_api.manager_service
import minknow_api.keystore_service

__all__ = [
    "Basecaller",
    "FlowCellPosition",
    "Manager",
]


class ServiceBase(object):
    def __init__(self, serviceclass, host, port, use_tls):
        self.host = host
        self.port = port
        if use_tls:
            self.channel = grpc.secure_channel(
                host + ":" + str(port), minknow_api.grpc_credentials()
            )
        else:
            self.channel = grpc.insecure_channel(host + ":" + str(port))
        self.rpc = serviceclass(self.channel)
        self.stub = self.rpc._stub

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if self.channel is not None:
            self.channel.close()
        self.rpc = None
        self.stub = None
        self.channel = None


class Manager(ServiceBase):
    """A connection to the manager gRPC interface.

    Args:
        host (str, optional): the hostname to connect to (note that IP addresses will not work for
            TLS connections)
        port (int, optional): override the port to connect to
        use_tls (bool, optional): set to False to use an insecure channel connection

    Attributes:
        bream_version (str): The version of Bream that is installed.
        config_version (str): The version of ont-configuration (Wanda) that is installed.
        channel (grpc.Channel): The gRPC channel used for communication.
        core_version (str): The running version of MinKNOW Core.
        core_version_components (tuple): A tuple of three integers describing the major, minor and
            patch parts of the core version. Useful for version comparisions.
        guppy_version (str): The version of Guppy that is running.
        host (str): The hostname used to connect.
        port (int): The port used to connect.
        rpc (minknow_api.manager_service.ManagerService): The auto-generated API wrapper.
        stub (minknow_api.manager_grpc_pb2.ManagerServiceStub): The gRPC-generated stub.
        version (str): The version of the MinKNOW distribution
        version_status (str): The status of the installed distribution ("unknown", "stable",
            "unstable" or "modified").
    """

    # The thread that handles some of these RPCs can block for a second or two in an operating
    # system call to start a process, especially on Windows. 5 seconds gives plenty of margin for
    # interference from antivirus hooks, for example, while still not blocking forever if something
    # has gone wrong with the manager.
    DEFAULT_TIMEOUT = 5

    def __init__(self, host="localhost", port=None, use_tls=True):
        if port is None:
            if use_tls:
                port = 9502
            else:
                port = 9501
        super(Manager, self).__init__(
            minknow_api.manager_service.ManagerService,
            host=host,
            port=port,
            use_tls=use_tls,
        )
        self._default_use_tls = use_tls

        version_info = self.stub.get_version_info(
            minknow_api.manager_service.GetVersionInfoRequest()
        )
        self.bream_version = version_info.protocols
        self.config_version = version_info.configuration
        self.core_version = version_info.minknow.full
        self.core_version_components = (
            version_info.minknow.major,
            version_info.minknow.minor,
            version_info.minknow.patch,
        )
        self.guppy_version = version_info.guppy_connected_version
        self.version = version_info.distribution_version
        DistributionStatus = (
            minknow_api.instance_pb2.GetVersionInfoResponse.DistributionStatus
        )
        self.version_status = DistributionStatus.Name(
            version_info.distribution_status
        ).lower()

    def __repr__(self):
        return "Manager({!r}, {!r})".format(self.host, self.port)

    def keystore(self):
        """
        Find the keystore service running for this manager.

        Returns:
            KeyStore: The gRPC service for the manager level keystore.
        """
        return minknow_api.keystore_service.KeyStoreService(self.channel)

    def basecaller(self, use_tls=None, timeout=DEFAULT_TIMEOUT):
        """Connect to the basecalling interface.

        Args:
            use_tls (bool, optional): Force use of TLS (True) or insecure (False) connection.
                By default, it will use the same type of channel as the manager connection, or
                whatever is available to connect on.
            timeout (float, optional): The maximum time to wait for the call to complete. Should
                usually be left at the default.

        Returns:
            minknow_api.manager.Basecaller: The wrapper for the Basecaller service, or None if the
                connection couldn't be made.
        """
        bc_api = self.stub.basecaller_api(
            minknow_api.manager_service.BasecallerApiRequest(), timeout=timeout
        )
        if use_tls is None:
            if bc_api.secure == 0 and bc_api.insecure != 0:
                use_tls = False
            elif bc_api.insecure == 0 and bc_api.secure != 0:
                use_tls = True
            else:
                use_tls = self._default_use_tls

        if use_tls:
            port = bc_api.secure
        else:
            port = bc_api.insecure
        if port == 0:
            return None
        return Basecaller(self.host, port, use_tls)

    def create_directory(self, name, parent_path="", timeout=DEFAULT_TIMEOUT):
        """Create a directory on the host.

        Args:
            name (str): The name of the directory to be created.
            parent_path (str, optional): The name of the directory to be created.
            timeout (float, optional): The maximum time to wait for the call to complete. Should
                usually be left at the default.

        Returns:
            str: the name of the created directory
        """
        request = minknow_api.manager_service.CreateDirectoryRequest(
            parent_path=parent_path, name=name
        )
        return self.stub.create_directory(request, timeout=timeout).path

    def describe_host(self, timeout=DEFAULT_TIMEOUT):
        """Get information about the machine running MinKNOW.

        Args:
            timeout (float, optional): The maximum time to wait for the call to complete. Should
                usually be left at the default.

        Returns:
            minknow_api.manager_service.DescribeHostResponse: The information about the host.
        """
        return self.stub.describe_host(
            minknow_api.manager_service.DescribeHostRequest(), timeout=timeout
        )

    def reset_position(self, position, force=False, timeout=DEFAULT_TIMEOUT):
        """Reset a flow cell position.

        Args:
            position (str): The name of the position to reset.
            force (bool, optional): Restart the position even if it seems to be running fine.
            timeout (float, optional): The maximum time to wait for the call to complete. Should
                usually be left at the default.
        """
        self.reset_positions([position], force=force, timeout=timeout)

    def reset_positions(self, positions, force=False, timeout=DEFAULT_TIMEOUT):
        """Reset flow-cell positions

        Args:
            positions (list of strings): The names of the positions to reset.
            force (bool, optional): Restart the position even if it seems to be running fine.
            timeout (float, optional): The maximum time to wait for the call to complete. Should
                usually be left at the default.
        """
        request = minknow_api.manager_service.ResetPositionRequest()
        request.positions.extend(positions)
        request.force = force
        self.stub.reset_position(request, timeout=timeout)

    def flow_cell_positions(self, timeout=DEFAULT_TIMEOUT):
        """Get a list of flow cell positions.

        Args:
            timeout (float, optional): The maximum time to wait for the call to complete. Should
                usually be left at the default.

        Yields:
            FlowCellPosition: A flow cell position. Ordering is not guaranteed.
        """
        call = self.stub.flow_cell_positions(
            minknow_api.manager_service.FlowCellPositionsRequest(), timeout=timeout
        )
        # avoid holding open the call for longer than necessary by consuming all the results up
        # front
        messages = [msg for msg in call]
        for msg in messages:
            for position in msg.positions:
                yield FlowCellPosition(position, self.host, self._default_use_tls)


class FlowCellPosition(object):
    """A flow cell position.

    Args:
        description (minknow_api.manager_service.FlowCellPosition): A description of a flow cell
            position returned from a call to ``flow_cell_positions`` or
            ``watch_flow_cell_positions`` on the manager.
        host (string, optional): The hostname of the manager API (see Manager.host). This will be
            used by the `connect` method. Defauls to "localhost".
        use_tls (bool, optional): Set the default for the ``use_tls`` argument for `connect`.
    """

    def __init__(self, description, host="localhost", use_tls=True):
        self.host = host
        self.description = description
        self._default_use_tls = use_tls
        self._device = None

    def __repr__(self):
        return "FlowCellPosition({!r}, {{{!r}}})".format(self.host, self.description)

    def __str__(self):
        return "{} ({})".format(self.name, self.state)

    @property
    def name(self):
        """str: The name of the position."""
        return self.description.name

    @property
    def location(self):
        """minknow_api.manager_service.FlowCellPosition.Location: The location of the position (if
        built-in, otherwise None).
        """
        if self.description.HasField("location"):
            return self.description.location
        else:
            return None

    @property
    def shared_hardware_group(self):
        """minknow_api.manager_service.FlowCellPosition.SharedHardwareGroup: The information about
        shared hardware (if built-in, otherwise None).
        """
        if self.description.HasField("shared_hardware_group"):
            return self.description.shared_hardware_group
        else:
            return None

    @property
    def state(self):
        """str: The state of the position.

        One of "initialising", "running", "resetting", "hardware_removed", "hardware_error" or
        "software_error".
        """
        State = minknow_api.manager_service.FlowCellPosition.State
        return State.Name(self.description.state)[6:].lower()

    @property
    def running(self):
        """bool: Whether the software for the position is running.

        Note that this is not directly equivalent to the "running" value of `state`: even when there
        are hardware errors, the software may still be running.
        """
        return self.description.HasField("rpc_ports")

    def connect(self, use_tls=None):
        """Connect to the position.

        Only valid to do if `running` is True.

        Returns:
            minknow_api.Connection: A connection to the RPC interface.
        """
        if not use_tls:
            use_tls = self._default_use_tls
        if use_tls:
            port = self.description.rpc_ports.secure
        else:
            port = self.description.rpc_ports.insecure
        return minknow_api.Connection(host=self.host, port=port, use_tls=use_tls)


class Basecaller(ServiceBase):
    """A connection to the basecalling gRPC interface.

    Args:
        host (str, optional): the hostname to connect to (note that IP addresses will not work for
            TLS connections)
        port (int, optional): override the port to connect to
        use_tls (bool, optional): set to False to use an insecure channel connection

    Attributes:
        channel (grpc.Channel): the gRPC channel used for communication
        host (str): the hostname used to connect
        port (int): the port used to connect
        rpc (minknow_api.manager_service.ManagerService): the auto-generated API wrapper
        stub (minknow_api.manager_grpc_pb2.ManagerServiceStub): the gRPC-generated stub
    """

    def __init__(self, host, port, use_tls=True):
        super(Basecaller, self).__init__(
            minknow_api.basecaller_service.Basecaller,
            host=host,
            port=port,
            use_tls=use_tls,
        )

    def __repr__(self):
        return "Basecaller({!r}, {!r})".format(self.host, self.port)
