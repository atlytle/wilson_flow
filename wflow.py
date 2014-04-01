import sys
from read_su3 import readfort, reconstruct

def main(files):
    dat = readfort(files[0])
    dat = reconstruct(dat) 
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
