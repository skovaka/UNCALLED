"""
A minimal server implementation that is compatible with
the minkow_api.Connection class. To use this server for
tests you should provide Servicer classes using keyword
arguments when initialising the server. Services should
implement all the methods that your tests require.

Eg:
>>> import minknow_api
>>> # Create a custom service
>>> class InstanceService(minknow_api.instance_pb2_grpc.InstanceServiceServicer):
...     def get_version_info(self, _request, _context):
...         pass
>>> # Add the custom service to the mock server
>>> server = minknow_api.testutils.MockMinKNOWServer(
...     instance_service=InstanceService
... )

There is one required method for the connection class
InstanceService.get_version_info() which supplies the
current version information. This is already included
in the mock server, however if you require methods in
the instance service then the mock InstanceService is
importable from this module:

>>> from minknow_api.testutils import MockMinKNOWServer, InstanceService
>>> class MyInstanceService(InstanceService):
...     def get_disk_space_info(self, request, context):
...         pass
>>> # Add the custom service to the mock server
>>> server = MockMinKNOWServer(instance_service=MyInstanceService)

For examples of other service implementations check the
test cases for this mock server.
"""
import logging
from concurrent import futures
from importlib import import_module
from packaging.version import parse

import grpc
import minknow_api
from minknow_api import _services, _optional_services

LOGGER = logging.getLogger(__name__)
VERSION = parse(minknow_api.__version__)
DEFAULT_SERVER_PORT = 0


class InstanceService(minknow_api.instance_pb2_grpc.InstanceServiceServicer):
    def get_version_info(self, _request, _context):
        """Find the version information for the instance"""
        return minknow_api.instance_pb2.GetVersionInfoResponse(
            minknow=minknow_api.instance_pb2.GetVersionInfoResponse.MinknowVersion(
                major=VERSION.major,
                minor=VERSION.minor,
                patch=VERSION.micro,
                full=minknow_api.__version__,
            )
        )


class MockMinKNOWServer:
    """A MinKNOW server that is compatible with the minknow_api.Connection

    This is a minimal gRPC server that implements all the required methods
    for the minknow_api.Connection class. The server is easily extensible
    using custom service classes and provides an interface for testing code
    that interacts with the MinKNOW without it being installed.

    Any documented service in the minknow_api module can be added to the
    server. The custom service should be passed as a keyword argument when
    initialising the server in the format ``{name}_service=MyClass`` where
    ``name`` is a valid minknow service.
    """

    def __init__(self, port=DEFAULT_SERVER_PORT, **kwargs):
        # Logging setup
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(name)s %(levelname)s %(message)s",
        )
        self.logger = logging.getLogger(__name__)

        # Set the server port or get an available port
        self.port = port
        if self.port == 0:
            self.logger.info("Port will be assigned by operating system")

        # Init the server
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))

        for service in kwargs:
            # get service name excluding '_service'
            if service.endswith("_service"):
                service = service[: -len("_service")]
            if service in _services:
                self.logger.info("Using user defined service for {!r}".format(service))
            else:
                self.logger.warning(
                    "Skipped user defined service {!r}, not recognised".format(service)
                )

        # Store kwargs, add 'instance_service' here as it is the
        #   only service/method that minknow.Connection requires
        self.kwargs = kwargs
        if "instance_service" not in self.kwargs:
            self.kwargs["instance_service"] = InstanceService

        # Iterate the services defined in minknow_api._services
        #   adding the base Servicer unless another is specified
        #   via kwargs
        for name, svc_list in _services.items():
            svc_name = "{}_service".format(name)
            try:
                svc_module = import_module("minknow_api.{}_pb2_grpc".format(name))
            except ImportError:
                if name not in _optional_services:
                    raise
            else:
                # There should be only one entry for each service
                for svc in svc_list:
                    svc_servicer = "{}Servicer".format(svc)

                    # Get the user overridden module from kwargs
                    #   or fallback onto the baseclass from grpc
                    module = self.kwargs.get(
                        svc_name, getattr(svc_module, svc_servicer),
                    )

                    # Same as: self.{name}_service = minknow_api.{name}_pb2_grpc.{svc}Servicer
                    #   or uses user defined module()
                    setattr(self, svc_name, module())

                    # Add servicer to server
                    func = "add_{}_to_server".format(svc_servicer)
                    add_servicer_to_server_func = getattr(svc_module, func)
                    add_servicer_to_server_func(getattr(self, svc_name), self.server)

        bound = self.server.add_insecure_port("[::]:{}".format(self.port))
        if bound == 0:
            raise ConnectionError("Could not connect using port {}".format(self.port))
        self.port = bound

        self.logger.info("Using port {}".format(self.port))

    def __enter__(self):
        self.start()
        return self.server

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop(0)

    def __getattr__(self, item):
        """Delegate attribute access to the gRPC server"""
        # Is this a bad idea?
        return getattr(self.server, item)
