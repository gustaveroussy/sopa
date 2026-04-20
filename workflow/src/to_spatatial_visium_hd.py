import argparse
import sopa.io

parser = argparse.ArgumentParser()

parser.add_argument("--input_path",
                     type = str,
                     help = 'path to the output of Space ranger',
                     required = True)

parser.add_argument("--output_path",
                     type = str,
                     help = 'path to the output of sopa to save .explorer and .zarr folders',
                     required = False)

parser.add_argument("--fullres_image",
                     type = str,
                     help = 'path to the fullres image file',
                     required = False)


sdata = sopa.io.visium_hd(path = args.input_path, fullres_image_file = args.fullres_image,dataset_id = "")

sdata.write(f"{args.output_path}", overwrite=True)