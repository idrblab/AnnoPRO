from annopro import main

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Arguments for main.py')
    parser.add_argument('--file_path', default=None, type=str)
    parser.add_argument('--used_gpu', default="0", type=str)
    parser.add_argument('--with_diamond', action='store_true', default=True)
    args = parser.parse_args()
    main(
        file_path=args.file_path,
        used_gpu=args.used_gpu,
        with_diamond=args.with_diamond
    )
