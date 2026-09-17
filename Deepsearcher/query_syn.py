import os
from dotenv import load_dotenv

load_dotenv()

from deepsearcher.configuration import Configuration, init_config
from deepsearcher.online_query import query

config = Configuration()

# Customize your config here,
# more configuration see the Configuration Details section below.
config.set_provider_config("llm", "DeepSeek", {"model": "deepseek-reasoner"})
config.set_provider_config("embedding", "SiliconflowEmbedding", {"model": "Pro/BAAI/bge-m3"})
init_config(config = config)

# Load your local data
from deepsearcher.offline_loading import load_from_local_files
load_from_local_files(paths_or_directory="/home/user/snow/deepsearch-z-25-3-24/ds/deep-searcher-master/paper")


# Query
result = query("How to design a synthesis protocol for a high-entropy Prussian blue analog (HE-PBA) catalyst using cobalt cyanide (KCo(CN)₆) as the framework and incorporating five transition metals (Fe, Co, Ni, Mn, Zn)? Please provide step-by-step instructions.",max_iter=10) # Your question here
