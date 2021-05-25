import uncalled as unc
import numpy as np

def main():
    parser = unc.ArgParser(
        "Rapidly maps raw nanopore signal to dna references", [
            unc.index.CMD, 
            unc.map.CMD, 
            unc.realtime.CMD, 
            unc.sim.CMD, 
            unc.pafstats.CMD
    ])

    cmd, conf = parser.parse_args()

    if cmd is not None:
        cmd(conf)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
