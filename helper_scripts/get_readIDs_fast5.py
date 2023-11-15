import sys
import glob
from ont_fast5_api.fast5_interface import get_fast5_file

fast5dir = sys.argv[1]


for file in glob.glob(fast5dir + '/*.fast5'):
    with get_fast5_file(file, mode="r") as f5:
        try:
            for read in f5.get_reads():
                print(read.read_id)
        except OSError as err:
            print("OS error: {0}".format(err))
            print("file: ", f5)
        except Exception as e:
            print(e)

# For reading a file
# with get_fast5_file(sys.argv[1], mode="r") as f5:
# 	print(f5)
# 	for read in f5.get_reads():
# 		print(read.read_id)

        #print("read ID:", read.read_id)
        #print("raw data", read.get_raw_data(scale=True))

        # ### only parse reads that are long enough
        # print("length: ", len(raw_data))
        # if len(raw_data) >= (cutoff + 3000):
        # 	print("if1")
        # 	print(read.read_id, "----------------------------")
        # 	if read.read_id in posli:
        # 		print("if2")
        # 		pi += 1
        # 		arrpos.append(raw_data[cutoff:(cutoff + 3000)])
        # 		print(pi, batch)
        # 		if (pi%batch == 0) and (pi != 0):
        # 			print("posi")
        # 			normalization(arrpos, pi, outpath, pos = True)
        # 			del arrpos
        # 			arrpos = []
