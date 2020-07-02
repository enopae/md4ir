# md4ir

A Python 3 code for treating MD data for the generation and analysis of IR spectra

## Basic usage

Call python on the `md4ir.py` code with `-h` flag to see the options. 
```
python md4ir.py -h
```

For each option, you can see the respective flags by invoking the option and calling the help feature again. For example, for spectrum calculation options:
```
python md4ir.py spec -h
```
### Usage recommendation

Generate a soft link to the python code in your `bin` directory:
```
ln -s path/to/md4ir.py md4ir
chmod 755 md4ir
```

You can now call `md4ir` in any directory as just:
```
md4ir -h
```


### Misc
The code has not been fully optimized. Disruptive updates are possible.

"Manual" to be updated soon...
