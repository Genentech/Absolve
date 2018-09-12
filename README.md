_____________________________________________________________________________
_____________________________________________________________________________

# Absolve

Absolve is a command line tool for annotating antibody variable domain sequences with 

	*  ORF
	*  Framework and CDR regions
	*  Kabat Numbering
	*  Germline Assignment
	*  Somatic Hypermutation 

Absolve accepts fasta files as input, or fastq files as either single or paired-end reads.  It generates results for the input sequences as either a single tab-delimited file or multiple tab-delimited files if requested.  Germline sequences are provided for human, but users maybe provide their own fasta formatted database files that have been indexed by BWA.  

Absolve is written in C++ and links directly to BWA and HMMR libraries included.

_____________________________________________________________________________
_____________________________________________________________________________

# Documentation

Please examine the manual for more detailed documentation of program parameters and design

## ./docs/manual/Absolve.html 

_____________________________________________________________________________
_____________________________________________________________________________

# Build the application

Supported platforms:
	* Linux
	* Mac OSX

1. Platform/OS dependencies:
	* C++11 compatible compiler
	* C mathematical library. (libm) 
	* Compression library. (libz)
	* POSIX multithreading library. (libpthread) 
	* Specific Boost libraries: boost_system, boost_filesystem, boost_iostreams, boost_regex
	* CMake Ver 3.5 or greater
	
2. Build the included dependencies.  Copies of BWA and HMMR are included.
```
	cd dep/
	tar -zxvf bwa-0.7.12.tar.gz
	tar -zxvf hmmer-3.1b2.tar.gz
	cd bwa-0.7.12/ 
	make
	cd ..
	cd hmmer-3.1b2/
	./configure
	make 
	cd ../..

```


3. Configure environment 
```
		source ./SETUP.bash `pwd`
```

4. Build release version 
```
		cmake CMakeLists.txt
		make clean
		make
```

5. Make sure all tests run ok 
```
		cd test
		cp -R testrun.ref.template testrun.ref
		./RunTest.bash
		cd ..
```

6.	Invoke *"absolve.bash"* to run

_____________________________________________________________________________
_____________________________________________________________________________
eof
