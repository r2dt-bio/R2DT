import select
import subprocess

from rich.console import Console
from rich import print


class Runner:
    def __init__(self, print_command: bool, print_output: bool):
        self.print_command = print_command
        self.print_output = print_output
        self.console = Console()

    def run(self, cmd: str) -> int:
        if self.print_command:
            print(f"[italic green]Executing:[/italic green] {cmd}")

        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

        while True:
            reads = [process.stdout.fileno(), process.stderr.fileno()]
            ret = select.select(reads, [], [])

            for fd in ret[0]:
                if fd == process.stdout.fileno():
                    line = process.stdout.readline()
                    if line:
                        self._stdout_callback(line.strip())
                if fd == process.stderr.fileno():
                    line = process.stderr.readline()
                    if line:
                        self._stderr_callback(line.strip())

            if process.poll() is not None:
                break


        # stdout, stderr = process.communicate()
        #
        # if self.print_output:
        #     if stdout:
        #         self.console.log(f"[italic green]STDOUT:[/italic green] {stdout.decode()}")
        #     if stderr:
        #         self.console.log(f"[italic red]STDERR:[/italic red] {stderr.decode()}")

        return process.returncode

    def _stdout_callback(self, line: str) -> None:
        if self.print_output:
            print(line)

    def _stderr_callback(self, line: str) -> None:
        if self.print_output:
            print(f"[red]{line}[/red]")


runner = Runner(print_command=True, print_output=True)
