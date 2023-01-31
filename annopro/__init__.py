
import sys
import os

if sys.platform == "win32":
    # numpy.distutils.command.build_ext.build_ext
    # will add dynamic library to the root package's .libs
    # because they think it's a good idea to share the same .libs
    # so we have to add it to the path
    # here is just for profeat fortran extension module
    shared_libs = os.path.join(os.path.dirname(__file__), ".libs")
    # During building, the .libs folder does not exist yet
    if os.path.isdir(shared_libs):
        os.add_dll_directory(shared_libs)


def console_main():
    import argparse
    parser = argparse.ArgumentParser(description='Arguments for AnnoPRO')
    parser.add_argument("--fasta_file", "-i", help="The protein sequences file")
    parser.add_argument('--output', "-o", default=None,
                        type=str, help="Output directory")
    parser.add_argument('--used_gpu', default="-1", type=str,
                        help="GPU device selected, default is CPU")
    parser.add_argument('--disable_diamond',
                        action='store_true', default=False,
                        help="Disable blast with diamond")
    parser.add_argument('--overwrite',
                        action="store_true",
                        default=False,
                        help="Overwrite existed output"
                        )
    parser.add_argument("--version",
                        action="store_true", default=False, help="Show version")
    args = parser.parse_args()
    if args.version:
        print("{} {}, Copyright Zhejiang University.".format(
            __name__, __version__))
        exit(0)
    elif args.fasta_file is None:
        parser.print_help()
        exit(1)
    main(
        proteins_fasta_file=args.fasta_file,
        output_dir=args.output,
        used_gpu=args.used_gpu,
        with_diamond=(not args.disable_diamond),
        overwrite=args.overwrite
    )


def main(proteins_fasta_file: str, output_dir: str = None,
         used_gpu: str = None, with_diamond: bool = True, overwrite: bool = False):
    from annopro.data_procession import profeat, process
    from diamond4py import Diamond
    from annopro import resources
    from os.path import join, exists
    from annopro.prediction import predict
    from shutil import rmtree

    if output_dir is None:
        output_dir = proteins_fasta_file + ".output"

    if exists(output_dir):
        if overwrite:
            rmtree(output_dir)
        else:
            print(f"Output directory {output_dir} already existed!")
            exit(1)

    profeat.run(proteins_fasta_file, output_dir)

    diamond_scores_file: str = None
    if with_diamond:
        diamond_scores_file = join(output_dir, "diamond_scores.txt")
        diamond = Diamond(
            database=resources.get_resource_path("cafa4.dmnd"),
            n_threads=4
        )
        diamond.blastp(
            query=proteins_fasta_file,
            out=diamond_scores_file
        )

    promap_features_file = join(output_dir, "promap_features.pkl")
    process(
        proteins_fasta_file=proteins_fasta_file,
        profeat_file=join(output_dir, "output-protein.dat"),
        save_file=promap_features_file)
    predict(output_dir=output_dir,
            promap_features_file=promap_features_file,
            used_gpu=used_gpu,
            diamond_scores_file=diamond_scores_file)


from . import _version
__version__ = _version.get_versions()['version']
