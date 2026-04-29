"""
Minimal stub for the original ``pyteltools.conf`` module.

Only the attributes consumed by the vendored `slf.Serafin` module are exposed
(LANG, ENABLE_MESH_ORIGIN). The upstream module relied on ``simple-settings``
to layer a default file with optional user overrides; here we replace that
machinery with a plain class so that no extra runtime dependency is needed.
"""


class _Settings:
    LANG = "fr"
    ENABLE_MESH_ORIGIN = True


settings = _Settings()
