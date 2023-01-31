import ipyrad as ip
import ipyrad.analysis as ipa
import ipyparallel as ipp
import sys


# if I decide to parallelize the runs
# Open a terminal and type the following command to start
# an ipcluster instance with 40 engines:
# ipcluster start -n 40 --cluster-id="ipyrad" --daemonize

# After the cluster is running you can attach to it with ipyparallel
# ipyclient = ipp.Client(cluster_id="ipyrad")

## the assembly's name
data = ip.Assembly(sys.argv[1])



# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()
