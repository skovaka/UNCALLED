import grpc

_wrapper_types = {
    "google.protobuf.BoolValue",
    "google.protobuf.BytesValue",
    "google.protobuf.DoubleValue",
    "google.protobuf.FloatValue",
    "google.protobuf.Int32Value",
    "google.protobuf.Int64Value",
    "google.protobuf.StringValue",
    "google.protobuf.UInt32Value",
    "google.protobuf.UInt64Value",
}


class MessageWrapper(object):
    def __init__(self, message, unwraps=[]):
        self._unwraps = unwraps
        self._objs = [getattr(message, attr) for attr in unwraps]
        self._objs.append(message)
        self._message = message

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return MessageWrapper(next(self._message))
        except grpc.RpcError as e:
            if e.code() == grpc.StatusCode.CANCELLED:
                raise StopIteration
            else:
                raise

    def __dir__(self):
        names = set([n for o in self._objs for n in dir(o)])
        names = names.union(dir(object))
        names = names.union(self.__dict__)
        return list(names)

    def __repr__(self):
        return "MessageWrapper({})".format(repr(self._message))

    def __str__(self):
        return str(self._message)

    def __eq__(self, other):
        if type(other) == MessageWrapper:
            return self._message == other._message
        else:
            return self._message == other

    def __ne__(self, other):
        if type(other) == MessageWrapper:
            return self._message != other._message
        else:
            return self._message != other

    def __iter__(self):
        for m in self._message:
            yield MessageWrapper(m, self._unwraps)

    def __getattr__(self, name):
        for o in self._objs:
            try:
                val = getattr(o, name)
            except AttributeError:
                # these aren't the droids you're looking for...
                pass
            else:
                try:
                    if val.DESCRIPTOR.full_name in _wrapper_types:
                        return val.value
                    else:
                        return val
                except AttributeError:
                    # clearly not a protobuf message
                    return val
        raise TypeError(
            '"{}" object has no attribute "{}"'.format(
                type(self._message).__name__, name
            )
        )


class ArgumentError(grpc.RpcError):
    """Exception raised by client-side wrapper code when argument processing fails."""

    def code(self):
        """Get the gRPC status code associated with this error."""
        return grpc.StatusCode.INVALID_ARGUMENT

    def details(self):
        """Get the gRPC error details associated with this error."""
        return str(self)
