import sys
from read_su3 import readfort, reconstruct
from su3 import plaq

def main(files):
    dat = readfort(files[0])
    dat = reconstruct(dat)
    print plaq(dat) 
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
