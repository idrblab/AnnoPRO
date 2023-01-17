import shutil
import os


def blastp(db: str, query: str, out: str):
    cmd_path = shutil.which("diamond")
    assert cmd_path is not None, "Please install diamond and add it to PATH"
    status = os.system(f"{cmd_path} blastp --db \"{db}\" --query \"{query}\" --out \"{out}\"")
    if status != 0:
        raise RuntimeError("An error occurred durring blasting")
