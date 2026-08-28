"""Tests for loading the run metadata absorption correction needs.

The absorption workspace is built from a run's metadata. Raw event NeXus is
read with ``LoadEventNexus``; a processed NeXus - which is what a simulated run
and any saved workspace looks like - has no ``/entry`` group and must fall back
to ``LoadNexus``. Only if neither works do we go looking on the archive.

Mantid's simpleapi wraps loader failures in ``RuntimeError``, so a fallback
that only catches ``ValueError`` never actually runs. It also infers
``OutputWorkspace`` from the assignment target at the call site, which a helper
function does not have - so the name has to be passed explicitly.
"""

import unittest

from total_scattering.file_handling import load as load_module


class FakeLoader:
    """Records its calls and fails with whatever it was told to fail with."""

    def __init__(self, error=None, result="workspace"):
        self.error = error
        self.result = result
        self.calls = []

    def __call__(self, filename, **kwargs):
        self.calls.append((filename, kwargs))
        if self.error is not None:
            raise self.error
        return self.result


class TestLoadAbsMetadata(unittest.TestCase):

    def test_raw_event_nexus_is_read_metadata_only(self):
        event = FakeLoader(result="event_ws")
        nexus = FakeLoader()

        result = load_module.load_abs_metadata(
            "NOM_144975", event, nexus, lambda name: "/archive/" + name)

        self.assertEqual(result, "event_ws")
        self.assertEqual(event.calls, [(
            "NOM_144975",
            {"MetaDataOnly": True, "OutputWorkspace": "abs_metadata_input"},
        )])
        self.assertEqual(nexus.calls, [])

    def test_output_workspace_is_named_explicitly(self):
        """simpleapi cannot infer it inside a helper, so we must pass it."""
        event = FakeLoader(error=RuntimeError("no /entry"))
        nexus = FakeLoader(result="processed_ws")

        load_module.load_abs_metadata(
            "sim.nxs", event, nexus, lambda name: name,
            output_workspace="my_ws")

        self.assertEqual(event.calls[0][1]["OutputWorkspace"], "my_ws")
        self.assertEqual(nexus.calls[0][1]["OutputWorkspace"], "my_ws")

    def test_processed_nexus_falls_back_to_load_nexus(self):
        """A saved or simulated workspace has no '/entry' group."""
        event = FakeLoader(
            error=RuntimeError(
                "LoadEventNexus-v1: supplied group /entry does not exist"))
        nexus = FakeLoader(result="processed_ws")

        result = load_module.load_abs_metadata(
            "sim.nxs", event, nexus, lambda name: "/archive/" + name)

        self.assertEqual(result, "processed_ws")
        self.assertEqual(nexus.calls, [(
            "sim.nxs", {"OutputWorkspace": "abs_metadata_input"})])

    def test_value_error_also_falls_back(self):
        event = FakeLoader(error=ValueError("not an event file"))
        nexus = FakeLoader(result="processed_ws")

        result = load_module.load_abs_metadata(
            "sim.nxs", event, nexus, lambda name: "/archive/" + name)

        self.assertEqual(result, "processed_ws")

    def test_falls_through_to_the_archive_path(self):
        """Neither loader handled the bare name, so resolve it on the archive."""
        event = FakeLoader(error=RuntimeError("no /entry"))
        nexus = FakeLoader(error=RuntimeError("not a processed file"))
        resolved = []

        def resolve(name):
            resolved.append(name)
            return "/SNS/NOM/IPTS-1/nexus/NOM_144975.nxs.h5"

        with self.assertRaises(RuntimeError):
            load_module.load_abs_metadata("NOM_144975", event, nexus, resolve)

        self.assertEqual(resolved, ["NOM_144975"])
        self.assertEqual(
            event.calls[-1],
            ("/SNS/NOM/IPTS-1/nexus/NOM_144975.nxs.h5",
             {"MetaDataOnly": True, "OutputWorkspace": "abs_metadata_input"}))


if __name__ == '__main__':
    unittest.main()  # pragma: no cover
