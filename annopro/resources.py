"""
This module manages all required resources for annopro.
"""
from io import FileIO, TextIOWrapper
import os
import wget
import hashlib

RESOURCE_DIR = os.path.join(os.path.expanduser("~"), ".annopro/data")
os.makedirs(RESOURCE_DIR, exist_ok=True)

RESOURCE_DICT = {
    "cafa_train.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/cafa_train.pkl",
        "md5sum": "07d3e4334c31c914efec3f52cae5e498"
    },
    "cafa4_del.csv": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/cafa4_del.csv",
        "md5sum": "d00b71439084cb19b7d3d0d4fbbaa819"
    },
    "cafa4.dmnd": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/cafa4.dmnd",
        "md5sum": "a2a6cba9af26dbe1911e14a306db2712"
    },
    "data_grid.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/data_grid.pkl",
        "md5sum": "fb2d2d86a4bc21c6e60fac996b3a90d3"
    },
    "go.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/go.pkl",
        "md5sum": "8d7a975d38a4af670b0370f4ea722a2a"
    },
    "go.txt": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/go.txt",
        "md5sum": "1dae308468fa00ae6d5796fd22c65044"
    },
    "row_asses.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/row_asses.pkl",
        "md5sum": "bf9bb1eda744a60c381d19b275ac6f33"
    },
    "terms_biological_process.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/terms_biological_process.pkl",
        "md5sum": "e79cf5e006432c19606b8a482cf7ddfa"
    },
    "terms_cellular_component.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/terms_cellular_component.pkl",
        "md5sum": "8cf58075eba65e2bb710566e4ef93f42"
    },
    "terms_molecular_function.pkl": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/data/terms_molecular_function.pkl",
        "md5sum": "002c3696fbca0061402a68129b35dcd4"
    },
    "bp.h5": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/model_param/bp.h5",
        "md5sum": "7e19158e5252a70ff831f5f583b1c2ed"
    },
    "cc.h5": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/model_param/cc.h5",
        "md5sum": "73876beec9370ff56b58878cf4446d2c"
    },
    "mf.h5": {
        "url": "https://promap.oss-cn-hangzhou.aliyuncs.com/model_param/mf.h5",
        "md5sum": "f4fb632f553afeb45571a29e46286bb8"
    }
}


def md5sum(file_path: str) -> str:
    with open(file_path, "rb") as f:
        md5 = hashlib.md5()
        while True:
            data = f.read(65536)
            if not data:
                break
            md5.update(data)
        return md5.hexdigest()


def md5check(file_path: str, expected: str):
    return md5sum(file_path).startswith(expected)


def download_resource(name: str, overwrite: bool = False) -> str:
    if name in RESOURCE_DICT:
        resource = RESOURCE_DICT[name]
        path_name = os.path.join(RESOURCE_DIR, name)
        if os.path.exists(path_name):
            if overwrite or not md5check(path_name, resource["md5sum"]):
                os.remove(path_name)
            else:
                return path_name
        print(f"Download {name}...")
        wget.download(
            url=resource["url"],
            out=path_name)
        print(f"\nValidate md5sum of {name}...")
        if not md5check(path_name, resource["md5sum"]):
            raise RuntimeError(f"{name} do not pass md5 validation, please visit https://github.com/idrblab/AnnoPRO for help")
        return path_name
    else:
        raise FileNotFoundError(f"Invalid resource name: {name}")


def get_resource_path(name: str) -> str:
    return download_resource(name)


def open_binary(name: str) -> FileIO:
    file_path = get_resource_path(name)
    return open(file_path, "rb")


def open_text(name: str) -> TextIOWrapper:
    file_path = get_resource_path(name)
    return open(file_path, "rt")
