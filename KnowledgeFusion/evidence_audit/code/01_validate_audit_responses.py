"""Validate structured LLM-auditor responses and bind them to their requests."""

from __future__ import annotations

import argparse

from audit_common import AuditError, add_common_validation_arguments, validate_responses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_validation_arguments(parser)
    args = parser.parse_args()
    try:
        rows = validate_responses(args.requests, args.responses, args.output)
    except AuditError as exc:
        parser.error(str(exc))
    print(f"Validated {len(rows)} response(s): {args.output}")


if __name__ == "__main__":
    main()
