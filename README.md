ycomp
=====

*ycomp* is a tool for working with Y-chromosome data from YFull and FTDNA.

Run `ycomp -h` for information on how to use the program.

Installation
------------

This package and its command-line interface can be installed easily using *pip*; simply run the following from the command line (within the directory of the repository).

```sh
pip3 install .
```

Usage
-----

Run `ycomp -h` to get help with using the tool, or `ycomp [command] -h` for help with specific commands or subcommands.

To get started, you probably want to do something like the following.

1.	Create a new directory in which to store the database files and analysis, and make it the working directory.

	```sh
	mkdir -p [path_to_new_directory]
	cd [path to new directory]
	```

2.	Download one or more (sub)trees from YFull for the haplogroups that you are interested in.

	```sh
	ycomp tree download-yfull --hg [haplogroup]
	```

	Use the haplogroup name that appears on the Yfull website.

3.	Download and import SNP and STR kit data from FTDNA for the haplogroups that you are interested in.

	```sh
	ycomp snp fetch-ftdna [group_name]
	ycomp str fetch-ftdna [group_name]
	```

	The group name may be the name as it appears on the FTDNA website (with spaces) or in the URL (with spaces replaced by hyphens).

4.	Import SNP and STR kit data for YFull kits that you wish to compare.

	```sh
	ycomp snp add-yfull ~/Downloads/SNP_for_YF00000_20210101.csv
	ycomp str add-yfull ~/Downloads/STR_for_YF00000_20210101.csv
	```

	These files can be downloaded from the YFull website, once you are logged in.

5.	Run analysis on the SNP and STR databases to compare a given kit against other kits.

	```sh
	ycomp snp analyze -k [kit_number]
	ycomp str analyze -k [kit_number]
	```

	For example, if you wish to run the analyses against the kit you imported in the previous step, replace `[kit_number]` with `YF00000`.

	This will generate CSV files containing the resluts of the analysis (in the working directory, by default).
