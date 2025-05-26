import os.path

import numpy as np
import pandas as pd
from gensim.models import word2vec
from cuspy import ConfigUtils


class KeywordsExtractor:

    def __init__(self, config):
        self._tag_token_pairs_fp = config.tag_token_pairs_fp
        self._words_vec_fp = config.words_vec_fp
        self._keywords = config.keywords
        self._keywords_synonyms_fp = config.keywords_synonyms_fp
        self._min_score = config.min_score

    def get_words(self):
        tag_token_pairs_df = pd.read_csv(self._tag_token_pairs_fp, sep='\t', encoding='utf-8')
        for _, row in tag_token_pairs_df.iterrows():
            tag_token_pairs = eval(row.tag_token_pairs)
            if len(tag_token_pairs) == 0:
                continue
            tags, tokens = zip(*tag_token_pairs)
            yield tokens

    def _load_word_vectors(self):
        with open(self._words_vec_fp, 'r', encoding='utf-8')as f:
            for i, line in enumerate(f.readlines()):
                if i == 0:
                    continue
                line = line.strip()
                if len(line) == 0:
                    continue
                ls = line.split(' ')
                token = ' '.join(ls[:-400])
                vec = ls[-400:]
                # vec = [float(v) for v in ls[-400:]]
                yield token.lower(), vec

    def _get_vec_for_words(self, words: [str]):
        word_to_vec = {}
        for token, vec in self._load_word_vectors():
            for word in words:
                if word == token:
                    word_to_vec[word] = np.array([float(v) for v in vec])
            if len(word_to_vec) == len(words):
                return word_to_vec

    def get_words_vectors(self):
        sentences = list(self.get_words())
        model = word2vec.Word2Vec(sentences, vector_size=400, window=5, min_count=1)
        model.wv.save_word2vec_format(self._words_vec_fp, binary=False)

    def get_synonyms(self):
        keyword_to_vec = self._get_vec_for_words(self._keywords)
        synonyms = {'word': [], 'cos': [], 'keyword': []}
        for word, vec in self._load_word_vectors():
            vec = np.array([float(v) for v in vec])
            for k, vec_k in keyword_to_vec.items():
                cos = vec.dot(vec_k) / (np.linalg.norm(vec) * np.linalg.norm(vec_k))
                if cos >= self._min_score:
                    synonyms['word'].append(word)
                    synonyms['cos'].append(cos)
                    synonyms['keyword'].append(k)
                    continue
        syn_df = pd.DataFrame(synonyms)
        syn_df = syn_df.sort_values(by='cos', ascending=False)
        syn_df.to_csv(self._keywords_synonyms_fp, sep='\t', encoding='utf-8', index=False)

    def process(self):
        self.get_words_vectors()
        self.get_synonyms()


if __name__ == "__main__":
    ke = KeywordsExtractor(ConfigUtils.load_config('../config.json').proj_config)
    ke.process()
