The program requires internet connection to access IMGTV database.
It also requires the latest version of python and the libraries distinguished by the command import "library_name".
To run the scrypt, download it, open cmd, change directory to where the scrypt is saved and run it on python "python assembler2.py".


>The program is specifically made to assemble Sanger sequences of antibody fragments (scFv).
> VECTOR RADIOBUTTONS
  The nucleotide sequence may be sequenced from plasmids used in the laboratory (pDNL6, pGEX, pHygro...) or
  it may be sequenced from any custom template, with custom primers which can be typed in.
>LINKER SEQUENCE RADIOBUTTONS
 scFv(s) coming from an in-vitro display library may contain a linker spacing VL and VH chains.
 The sequence of the linker is utilized to distinguish sense and antisense strand Sanger sequences.
 This step is mandatory as sense and antisense strands are then merged based on nucleotide quality scores and further processed.
>hit ESEGUI! to run the assembler.
>After hitting the run button, the program requires .avi files from a directory. Multiple files can be selected. The program sorts sense and antisense sequences automatically.
>The program will generate a word file, where the output will be displayed.
>THE OUTPUT:
>The program will align forward and reverse sequence of each scFv, and merge them based on each the highest nucleotide base-call quality. 
 This will generate a complete forward sequence with the best quality.
>The alignment of forward and reverse strand will be displayed.
>The merged sequence will be aligned with human immunoglobulin genes, interrogating the IMGTV databse. The highest homology genes will be displayed in a table format
