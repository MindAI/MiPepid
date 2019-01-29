import numpy as np
import pickle
import pandas as pd

class kmer_featurization:

  def __init__(self, k):
    """
    k: the "k" in k-mer
    """
    self.k = k
    self.letters = ['A', 'T', 'C', 'G']
    self.multiplyBy = 4 ** np.arange(k-1, -1, -1) # the multiplying number for each digit position in the k-number system
    self.n = 4**k # number of possible k-mers

  def obtain_kmer_feature_for_a_list_of_sequences(self, seqs, write_number_of_occurrences=False):
    """
    Given a list of m DNA sequences, return a 2-d array with shape (m, 4**k) for the 1-hot representation of the kmer features.

    Args:
      write_number_of_occurrences:
        a boolean. If False, then in the 1-hot representation, the percentage of the occurrence of a kmer will be recorded; otherwise the number of occurrences will be recorded. Default False.    
    """
    kmer_features = []
    for seq in seqs:
      this_kmer_feature = self.obtain_kmer_feature_for_one_sequence(seq.upper(), write_number_of_occurrences=write_number_of_occurrences)
      kmer_features.append(this_kmer_feature)

    kmer_features = np.array(kmer_features)

    return kmer_features

  def obtain_kmer_feature_for_one_sequence(self, seq, write_number_of_occurrences=False):
    """
    Given a DNA sequence, return the 1-hot representation of its kmer feature.

    Args:
      seq: 
        a string, a DNA sequence
      write_number_of_occurrences:
        a boolean. If False, then in the 1-hot representation, the percentage of the occurrence of a kmer will be recorded; otherwise the number of occurrences will be recorded. Default False.
    """
    number_of_kmers = len(seq) - self.k + 1

    kmer_feature = np.zeros(self.n)

    for i in range(number_of_kmers):
      this_kmer = seq[i:(i+self.k)]
      this_numbering = self.kmer_numbering_for_one_kmer(this_kmer)
      kmer_feature[this_numbering] += 1

    if not write_number_of_occurrences:
      kmer_feature = kmer_feature / number_of_kmers

    return kmer_feature

  def kmer_numbering_for_one_kmer(self, kmer):
    """
    Given a k-mer, return its numbering (the 0-based position in 1-hot representation)
    """
    digits = []
    for letter in kmer:
      digits.append(self.letters.index(letter))

    digits = np.array(digits)

    numbering = (digits * self.multiplyBy).sum()

    return numbering

def predict(logr, X, threshold):
  y_pred_score = logr.predict_proba(X)[:,1]
  y_pred = (y_pred_score > threshold) + 0
  y_prob = abs(1 - y_pred - y_pred_score)  # the probability that an instance is in the assigned category
  return y_pred, y_pred_score, y_prob

def load_model(model_fname = './src/model/model.pkl'):
  f = open(model_fname, 'rb')
  logr = pickle.load(f)
  threshold = pickle.load(f)
  f.close()
  return logr, threshold

def predict_on_one_batch_and_write(sORFs, logr, threshold, output_fname, k=4):
  class_dic = {
    1: 'coding',
    0: 'noncoding'
  }

  columns=['sORF_ID', 'sORF_seq', 'transcript_DNA_sequence_ID', 'start_at', 'end_at'] # 'start_at' and 'end_at' are both 1-based.
  df = pd.DataFrame(sORFs, columns=columns)
  seqs = df.sORF_seq.tolist()
  obj = kmer_featurization(k)
  kmer_features = obj.obtain_kmer_feature_for_a_list_of_sequences(seqs, write_number_of_occurrences=False)
  y_pred, _, y_prob = predict(logr, kmer_features, threshold)
  df['classification'] = [class_dic[x] for x in y_pred]
  df['probability'] = y_prob
  df.to_csv(output_fname, header=False, index=False, mode='a')