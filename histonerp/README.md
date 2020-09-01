### Prerequisites

#### Mac

``` bash
brew install openssl
export C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/local/Cellar/openssl/your_version/include
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/Cellar/openssl/your_version/lib/"
```

#### Linux
``` bash
sudo apt-get install openssl
```

### Install the module using:

``` bash
sudo python setup.py install
# or 
python setup.py install --user
```

### Module on the production:

#### test the module:

``` bash
./getbw
``` 

``` bash
import _bw
_bw.getrp('your_bigWig', 'annotation_TSS', 'output', 10000.0)
```

#### Command line
regpot.py  -n out -b test/42359_treat.bw --tss test/test_tss.bed

adjust `-a ALPHA` for decay rate.

#### Further instruction
annotation_TSS `--tss` should be formatted as:

    chr1    11873   11874   NR_046018:DDX11L1       0       +
    chr1    17436   17437   NR_107063:MIR6859-3     0       -
    chr1    17436   17437   NR_128720:MIR6859-4     0       -

output demo `-n`, `5th` column is the regulatory potential score from bigWig: 

    chr1    11873   11874   NR_046018       208.234170501   DDX11L1 +
    chr1    17436   17437   NR_107063       219.513716641   MIR6859-3       -
    chr1    17436   17437   NR_128720       219.513716641   MIR6859-4       -


#### parameter to delete one specific region

This removal of region does not consider strand, e.g., to delete all promoter +/- 1kb:

``` bash
_bw.getrp('test/test.bw', 'test/test_tss.bed', 'test/out.tab', 10000.0, -1000, 1000)
```

Below would delete the downstream 1kb in the gene body for positive strand transcript, but upstream 1kb for negative strand transcript.

``` bash
_bw.getrp('test/test.bw', 'test/test_tss.bed', 'test/out.tab', 10000.0, 0, 1000)
```

### Compatible with both Python 2 and 3...

Remember the difference is so much between 2 and 3, be careful for: 

``` bash
1. PyInt is replaced with PyLong
2. PyString is replaced with PyBytes
3. Py_InitModule is replaced with PyModule_Create with PyModuleDef 
4. void initmodule(void) is replaced with PyMODINIT_FUNC PyInit_module(void)
Trick: use #if PY_MAJOR_VERSION >= 3 to define aliases for python3 
such as #define PyInt_Type PyLong_Type
``` 


Follow the demo here https://docs.python.org/release/3.0.1/howto/cporting.html#cporting-howto.
reference:
```
1. https://docs.python.org/3/extending/extending.html
2. http://python3-cookbook.readthedocs.io/zh_CN/latest/c15/p02_write_simple_c_extension_module.html
3. https://docs.python.org/3/c-api/long.html#c.PyLong_FromLong
```

