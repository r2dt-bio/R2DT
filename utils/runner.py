import select
import subprocess
from typing import Optional

from rich import print


class Runner:
    def __init__(self, print_command: bool, print_output: bool):
        self.print_command = print_command
        self.print_output = print_output

    def run(self, cmd: str, print_output: Optional[bool] = None) -> int:
        if self.print_command:
            print(f"[green]Executing:[/green] {cmd}")

        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

        while True:
            reads = [process.stdout.fileno(), process.stderr.fileno()]
            ret = select.select(reads, [], [])

            for fd in ret[0]:
                if fd == process.stdout.fileno():
                    line = process.stdout.readline()
                    if line and line.strip():
                        self._stdout_callback(line.strip(), print_output)
                if fd == process.stderr.fileno():
                    line = process.stderr.readline()
                    if line and line.strip():
                        self._stderr_callback(line.strip())

            if process.poll() is not None:
                break

        return process.returncode

    def _stdout_callback(self, line: str, print_output: Optional[bool]) -> None:
        if print_output is None and self.print_output:
            print(line, flush=True)
        elif print_output:
            print(line, flush=True)

    def _stderr_callback(self, line: str) -> None:
        if self.print_output or True:
            print(f"[red]E: {line}[/red]", flush=True)


runner = Runner(print_command=False, print_output=False)
