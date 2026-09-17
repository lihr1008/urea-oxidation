"""Build frozen KnowledgeFusion inputs from validated audit and adjudication data."""

from __future__ import annotations

import argparse
import json

from audit_common import AuditError, build_frozen_inputs, path_arg, verify_reference_frozen_inputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atomization", type=path_arg)
    parser.add_argument("--mapping", type=path_arg)
    parser.add_argument("--classification", type=path_arg)
    parser.add_argument("--adjudication", type=path_arg)
    parser.add_argument("--output-dir", type=path_arg)
    parser.add_argument(
        "--verify-reference",
        type=path_arg,
        help="Read-only validation of an existing frozen-input directory.",
    )
    args = parser.parse_args()
    try:
        if args.verify_reference:
            result = verify_reference_frozen_inputs(args.verify_reference)
        else:
            missing = [
                name
                for name, value in [
                    ("--atomization", args.atomization),
                    ("--mapping", args.mapping),
                    ("--classification", args.classification),
                    ("--adjudication", args.adjudication),
                    ("--output-dir", args.output_dir),
                ]
                if value is None
            ]
            if missing:
                raise AuditError(f"missing required arguments: {', '.join(missing)}")
            result = build_frozen_inputs(
                args.atomization,
                args.mapping,
                args.classification,
                args.adjudication,
                args.output_dir,
            )
    except AuditError as exc:
        parser.error(str(exc))
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
