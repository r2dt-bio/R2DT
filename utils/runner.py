import select
import subprocess
from typing import Optional

from rich import print as rprint


class Runner:
    """Executor for subprocesses with logging"""

    def __init__(self, print_command: bool, print_output: bool):
        self.print_command = print_command
        self.print_output = print_output
        self.always_print_stderr = True

    def run(self, cmd: str, print_output: Optional[bool] = None) -> int:
        """Execute a command in a subprocess and print results if needed"""

        if self.print_command:
            rprint(f"[green]Executing:[/green] {cmd}")

        with subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True,
        ) as process:
            while True:
                reads = [process.stdout.fileno(), process.stderr.fileno()]
                ret = select.select(reads, [], [])

                for file_descriptor in ret[0]:
                    if file_descriptor == process.stdout.fileno():
                        line = process.stdout.readline()
                        if line and line.strip():
                            self._stdout_callback(line.strip(), print_output)
                    if file_descriptor == process.stderr.fileno():
                        line = process.stderr.readline()
                        if line and line.strip():
                            self._stderr_callback(line.strip())

                if process.poll() is not None:
                    break

            return process.returncode

    def _stdout_callback(self, line: str, print_output: Optional[bool]) -> None:
        if print_output is None and self.print_output:
            rprint(line, flush=True)
        elif print_output:
            rprint(line, flush=True)

    def _stderr_callback(self, line: str) -> None:
        if self.print_output or self.always_print_stderr:
            rprint(f"[red]E: {line}[/red]", flush=True)


runner = Runner(print_command=False, print_output=False)
