class AnyMeta(type):
    def __getattr__(cls, k):
        return AnyMeta(cls.__name__ + "_" + k, (Any,), {})
class Any(metaclass=AnyMeta):
    def __init__(self, *a, **k): pass
    def __getattr__(self, k): return Any()
    def __call__(self, *a, **k): return Any()
