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
load_from_local_files(paths_or_directory="/home/user/snow/deepsearch-z-25-3-24/ds/deep-searcher-master/data")


# Query
result = query("In an alkaline environment for the urea oxidation reaction, design a high-entropy Prussian blue analog catalyst composed of five transition metals. This catalyst should not only exhibit high catalytic activity and stability but also meet multiple criteria including low cost, ease of synthesis, environmental friendliness, and corrosion resistance. Among common metal combinations, considering the balance between catalytic performance and structural stability, and aiming to avoid certain metals that may pose higher costs or environmental concerns, a tentative introduction of a metal element with stable valence states and environmentally benign characteristics is proposed to optimize catalyst performance. Please analyze in detail the rationale behind the selection of each of the five transition metal elements, elaborating on their contributions in terms of catalytic activity, stability, synthetic feasibility, cost-effectiveness, and structural features.",max_iter=8)
