"""Prepare client-neutral EvidenceAudit request batches without calling an LLM."""

from __future__ import annotations

import argparse

from audit_common import AuditError, path_arg, prepare_requests


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=["atomization", "mapping", "classification", "recheck"],
        required=True,
    )
    parser.add_argument("--run-packets", type=path_arg, required=True)
    parser.add_argument("--responses", type=path_arg)
    parser.add_argument("--revised-claims", type=path_arg)
    parser.add_argument("--output", type=path_arg, required=True)
    args = parser.parse_args()
    try:
        rows = prepare_requests(
            phase=args.phase,
            run_packets_path=args.run_packets,
            responses_path=args.responses,
            revised_claims_path=args.revised_claims,
            output_path=args.output,
        )
    except AuditError as exc:
        parser.error(str(exc))
    print(f"Prepared {len(rows)} {args.phase} request(s): {args.output}")


if __name__ == "__main__":
    main()
