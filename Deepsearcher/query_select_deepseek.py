import sys
import io

# Force the encoding of standard output and standard error to UTF-8 (Fixes UnicodeEncodeError)
# ⚠️ Note: This code MUST be placed at the very beginning of the file to ensure execution before any print statements.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

import os
from dotenv import load_dotenv

load_dotenv()

from deepsearcher.configuration import Configuration, init_config
from deepsearcher.online_query import query
from deepsearcher.offline_loading import load_from_local_files

# Import the new pymilvus client
from pymilvus import MilvusClient

COLLECTION_NAME = "deepsearcher"
MILVUS_URI = "http://localhost:19530"


# First, check Milvus status (using the more modern MilvusClient approach)
def check_and_load_data():
    # Connect to Milvus to check status
    try:
        # Use MilvusClient for connection and checking, maintaining consistency with Deepsearcher configuration
        client = MilvusClient(uri=MILVUS_URI)

        if client.has_collection(COLLECTION_NAME):
            print(f"✅ Collection '{COLLECTION_NAME}' already exists, skipping data loading")
            return False  # Data loading is not required
        else:
            print(f"🔄 Collection '{COLLECTION_NAME}' does not exist, data loading is required")
            return True  # Data loading is required
    except Exception as e:
        # Catch errors like connection failure
        print(f"❌ Milvus connection or check failed: {e}")
        # If connection fails, default to True, allowing load_from_local_files to attempt creating the connection and collection
        return True


config = Configuration()

config.set_provider_config("llm", "DeepSeek", {"model": "deepseek-reasoner"})
config.set_provider_config("embedding", "SiliconflowEmbedding", {"model": "Pro/BAAI/bge-m3"})

DATA_DIR = r"D:\pythonProject\program\urea\high entropy materials for urea oxidation\Deepsearch\data"
config.set_provider_config("vector_db", "Milvus", {
    "uri": MILVUS_URI,
    "collection_name": COLLECTION_NAME,
})

# Initialize configuration
init_config(config=config)

# Decide whether to load data based on the check result
need_load = check_and_load_data()

if need_load:
    print("Starting data loading...")
    try:
        load_from_local_files(paths_or_directory=DATA_DIR)
        print("✅ Data loading complete")
    except Exception as e:
        print(f"❌ Data loading failed: {e}")
else:
    print("Using existing data")


import nni
def get_default_parameters():
    # Default parameters for the decision tree model
    params = {
        'try_count': 12,  # 3,6
    }
    return params


# received_PARAMS = nni.get_next_parameter()
config2 = get_default_parameters()
# config2.update(received_PARAMS)

result = query(
    "In an alkaline environment for the urea oxidation reaction, design a high-entropy Prussian blue analog catalyst composed of five transition metals. This catalyst should not only exhibit high catalytic activity and stability but also meet multiple criteria including low cost, ease of synthesis, environmental friendliness, and corrosion resistance. Among common metal combinations, considering the balance between catalytic performance and structural stability, and aiming to avoid certain metals that may pose higher costs or environmental concerns, a tentative introduction of a metal element with stable valence states and environmentally benign characteristics is proposed to optimize catalyst performance. Please analyze in detail the rationale behind the selection of each of the five transition metal elements, elaborating on their contributions in terms of catalytic activity, stability, synthetic feasibility, cost-effectiveness, and structural features.",
    max_iter=30
)

# 将结果输出到txt文件
with open(f"deepseek_result/{config2['try_count']}_30.txt", "w", encoding="utf-8") as f:
    f.write(str(result))

print("查询结果已保存到 query_result.txt 文件中")