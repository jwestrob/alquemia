#!/bin/bash
# Login-host Claude API client/orchestrator. Do not sbatch this script.
# Arguments: WORKTREE TASK_SPEC OUTPUT_DIRECTORY MODEL EFFORT SESSION_UUID CANONICAL_REPO
set -euo pipefail
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo 'Claude sessions run on the login host, not inside Slurm.' >&2
  exit 2
fi
if [[ $# -ne 7 ]]; then
  echo 'Expected WORKTREE TASK_SPEC OUTPUT_DIRECTORY MODEL EFFORT SESSION_UUID CANONICAL_REPO' >&2
  exit 2
fi
TASK_WORKTREE="$1"
TASK_SPEC="$2"
TASK_OUTPUT="$3"
TASK_MODEL="$4"
TASK_EFFORT="$5"
TASK_SESSION_UUID="$6"
CANONICAL_REPO="$7"
SPEC_DIR="$(dirname "$TASK_SPEC")"
CLAUDE_BIN=/home/jwestrob/.local/bin/claude
for task_file in "$TASK_SPEC" "$SPEC_DIR/COMMON_RULES.md" "$SPEC_DIR/agents.json" "$SPEC_DIR/session.settings.json"; do
  test -f "$task_file"
done
test -d "$TASK_WORKTREE"
mkdir "$TASK_OUTPUT"
cd "$TASK_WORKTREE"
"$CLAUDE_BIN" --version > "$TASK_OUTPUT/cli_version.txt"
set +e
"$CLAUDE_BIN" auth status > "$TASK_OUTPUT/auth_status.json" 2> "$TASK_OUTPUT/auth.stderr"
auth_status=$?
set -e
if [[ "$auth_status" -ne 0 ]]; then
  echo 'blocked_auth' > "$TASK_OUTPUT/BLOCKED.txt"
  exit "$auth_status"
fi
{
  cat "$SPEC_DIR/COMMON_RULES.md"
  printf '\nCanonical repository (read-only artifacts): %s\nOwned worktree: %s\nOwned output: %s\n' "$CANONICAL_REPO" "$TASK_WORKTREE" "$TASK_OUTPUT"
  printf '\nBefore task work, read applicable CLAUDE.md AND AGENTS.md as specified in COMMON_RULES.md, then docs/AGENT_PIPELINE.md and the current field-response checkpoint. Record instructions_read paths in STATUS.json. Newer receipts supersede old running-status prose.\n'
  cat "$TASK_SPEC"
} > "$TASK_OUTPUT/prompt.txt"
sha256sum "$TASK_SPEC" "$SPEC_DIR/COMMON_RULES.md" "$SPEC_DIR/agents.json" "$SPEC_DIR/session.settings.json" > "$TASK_OUTPUT/spec_hashes.txt"
set +e
"$CLAUDE_BIN" -p --model "$TASK_MODEL" --effort "$TASK_EFFORT" \
  --session-id "$TASK_SESSION_UUID" --output-format stream-json --verbose \
  --permission-mode dontAsk --permission-prompts none \
  --settings "$SPEC_DIR/session.settings.json" --agents "$SPEC_DIR/agents.json" \
  --tools 'Read,Glob,Grep,Write,Edit,Bash,Agent' \
  --strict-mcp-config --mcp-config '{"mcpServers":{}}' \
  --disable-slash-commands --add-dir "$CANONICAL_REPO" "$TASK_OUTPUT" \
  < "$TASK_OUTPUT/prompt.txt" > "$TASK_OUTPUT/events.jsonl" 2> "$TASK_OUTPUT/stderr.log"
claude_status=$?
set -e
python3 - "$TASK_OUTPUT" "$TASK_SESSION_UUID" "$claude_status" <<'PY'
import datetime,json,os,pathlib,sys
out,session,code=sys.argv[1:]
p=pathlib.Path(out)/'session_exit.json'
p.write_text(json.dumps({'session_uuid':session,'launcher_pid':os.getppid(), 'host':os.uname().nodename,'returncode':int(code),'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task_accepted':False,'note':'Controller must inspect result event, tool denials, task artifacts and required tests.'},indent=2)+'\n')
PY
exit "$claude_status"
