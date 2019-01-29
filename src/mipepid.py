from ORF import ORFs, collect_and_name_sORFs_from_an_ORFs_object
import ML

from Bio import SeqIO
import pandas as pd
import sys

def MiPepid(input_fname, output_fname):
  # initialize the output file and write the header
  columns=['sORF_ID', 'sORF_seq', 'transcript_DNA_sequence_ID', 'start_at', 'end_at', 'classification', 'probability']
  df = pd.DataFrame(columns=columns)
  df.to_csv(output_fname, index=False)
  print('Begin writing the output file: ' + output_fname)

  # load the model
  logr, threshold = ML.load_model()

  # gather all the sORFs and write
  all_sORFs = []

  for rec in SeqIO.parse(input_fname, 'fasta'):
    DNA_seq = str(rec.seq).upper()
    obj_ORFs = ORFs(DNA_seq)
    transcript_seq_ID = rec.id
    this_sORFs = collect_and_name_sORFs_from_an_ORFs_object(obj_ORFs, transcript_seq_ID)
    all_sORFs += this_sORFs

    # Process in batch, 1 batch slightly greater than 1000
    if len(all_sORFs) > 1000:
      ML.predict_on_one_batch_and_write(all_sORFs, logr, threshold, output_fname)
      all_sORFs = []
      print('Wrote another 1000 sORFs.')

  if len(all_sORFs) > 0:
    ML.predict_on_one_batch_and_write(all_sORFs, logr, threshold, output_fname)
    print('Finished writing all the sORFs.')


if __name__=='__main__':
  input_fname = output_fname = None

  if len(sys.argv) == 1:
    print('Please specifiy the input fasta file using the following command:')
    print('python3 ./src/mipepid.py input_file_name')
  elif len(sys.argv) == 2:
    input_fname = sys.argv[1]
    output_fname = 'Mipepid_results.csv'
  else:
    input_fname = sys.argv[1]
    output_fname = sys.argv[2]

  if input_fname and output_fname:
    MiPepid(input_fname, output_fname)


  