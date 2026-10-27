from ._any import Any, AnyMeta
def __getattr__(name): return AnyMeta("core_" + name, (Any,), {})
