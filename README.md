# PILFER

PILFER (PIrna cLuster FindER) is a tool to predict clusters in piRNA sequences given their count and the loci. This tool reports clusters based on the assumption that the piRNAs are conserved and highly expressed in clusters. Also included the piRNA analysis pipeline used for the publication. You need to change the variables used in the script to make it work for your own data. Also, you need to build your own piRNA indexes n case you would like to use it on other species.

## Getting Started

The Python script pilfer.py is included in the zip archive, along with a sample input file in the BED format. Remember that the count value of each read should be present in the fifth column of the BED file as an integer value.

### Prerequisites

You need to have Python 2.7 in your system and the numpy module installed. Most Python distribution include the numpy module natively in the installation. If you are using Python 2.6 or below, please make sure the following packages are installed

* Numpy
* Operator
* Sys
* CSV
* OS
* Argparse

The script was tested under Linux environment.

## Installing and Running the script

All you need is the script to predict the clusters in your input file!

### Run the script

You can run the script with the -h or --help option to view the help

```
python pilfer.py --help
```
To run the script type

```
python pilfer.py -i sample_input.bed
```

To set a standard deviation factor, use the -f switch

```
python pilfer.py -f 3.0 -i sample_input.bed
```
###Output
By default the output is printed to the standard output. If you want to print it to a file, you can simply redirect it to a file using redirection operator '>'

```
python pilfer.py -i sample_input.bed > output.txt
```

## Built With

* Python 2.7

## Authors

* **Rishav Ray** - *Initial work* - [Github](https://github.com/rishavray) rishav.rray@gmail.com
With contribution of the mouse dataset and the pipeline from Giovanni Pascarella, Division of Genomic Medicine, RIKEN IMS, Japan. giovanni.pascarella@riken.jp

## Citation
Ray, Rishav, and Priyanka Pandey. "piRNA analysis framework from small RNA-Seq data by a novel cluster prediction tool-PILFER." Genomics (2017).

## License

This project is licensed under the MIT License - see the [LICENSE](https://opensource.org/licenses/MIT) file for details
