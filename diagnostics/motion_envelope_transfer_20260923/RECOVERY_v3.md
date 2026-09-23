# Read-only scalar cell join recovery

The successful finite candidate calculations were missing only the collector's
`cell_id` association. Their scalar tasks already contain exact case, metal and
candidate identities. The new `motion_envelope_transfer_collect.py` maps those keys
to one unique original pool task, checks charge/multiplicity/coordinates, and applies
the existing collector with that in-memory association. It writes separate JOIN
receipts and never modifies the original manifests, outputs, energy values or
execution receipts. The actual SCF failure remains unavailable.

Both scheduled pool collectors and their dependent report job failed before this
repair. All allocation costs remain included. Read-only recovery produced the
originally requested final collection filenames only after confirming those files
did not exist; no successful result was overwritten. The final comparator runs on
those collections with its already frozen implementation and reference. No molecular
calculation was rerun or added by this repair.
