"""
Helpers for accessing the device service
========================================

This provides a `DeviceType` enum, which is a more Pythonic way of dealing with the ``device_type``
field from ``device.get_device_info`` RPC. The `get_device_type` function provides a convenient way
to obtain this value.

"""

import minknow_api.device_service
from enum import Enum

__all__ = ["DeviceType", "get_device_type"]


class DeviceType(Enum):
    """The type of device."""

    MINION = minknow_api.device_service.GetDeviceInfoResponse.MINION
    PROMETHION = minknow_api.device_service.GetDeviceInfoResponse.PROMETHION
    GRIDION = minknow_api.device_service.GetDeviceInfoResponse.GRIDION
    MINION_MK1C = minknow_api.device_service.GetDeviceInfoResponse.MINION_MK1C
    TRAXION = minknow_api.device_service.GetDeviceInfoResponse.TRAXION

    def is_minion_like(self):
        """Whether the device acts like a MinION.

        Among other things, this means it can be used with the ``minion`` RPC service (see
        `minknow_api.minion_service`).
        """
        return self in [
            DeviceType.MINION,
            DeviceType.GRIDION,
            DeviceType.MINION_MK1C,
            DeviceType.TRAXION,
        ]

    def is_promethion_like(self):
        """Whether the device acts like a PromethION.

        Among other things, this means it can be used with the ``promethion`` RPC service (see
        `minknow_api.prometion_service`).
        """
        return self in [DeviceType.PROMETHION]


def get_device_type(connection):
    """Get a `DeviceType` enum value from a connection to a flow cell position.

    Args:
        connection (minknow_api.Connection): Connection to a MinKNOW flow cell position.

    Results:
        DeviceType: The type of flow cell position.
    """
    return DeviceType(connection.device.get_device_info().device_type)
