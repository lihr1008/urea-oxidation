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

# Please input your API KEY
API_KEY = os.getenv("API_KEY", "")
# 必须指定完整的 API URL 根目录
# CUSTOM_BASE_URL = "https://api.yyds168.net/v1/"
#CUSTOM_BASE_URL = "https://api.probex.top/v1/"
CUSTOM_BASE_URL="https://api.vectorengine.ai/v1/"

config.set_provider_config("llm", "OpenAI", {
    "model": "gpt-5",
    "base_url": CUSTOM_BASE_URL,  # 👈 修正：填入您的自定义 API 根 URL
    "api_key": API_KEY      # 👈 使用您有效的密钥
})

config.set_provider_config("embedding", "SiliconflowEmbedding", {"model": "Pro/BAAI/bge-m3"})

DATA_DIR = r"D:\pythonProject\program\urea\high entropy materials for urea oxidation\Deepsearch\paper"
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
        'try_count': 5,  # 25,57,93,16,42,131,197,107
    }
    return params


# received_PARAMS = nni.get_next_parameter()
config2 = get_default_parameters()
# config2.update(received_PARAMS)

result = query(
    "How to design a synthesis protocol for a high-entropy Prussian blue analog (HE-PBA) catalyst using cobalt cyanide K3[Co(CN)6] as the framework and incorporating five transition metals (Fe, Co, Ni, Mn, Zn)? Please provide step-by-step instructions.",
    max_iter=50
)

# 将结果输出到txt文件
with open(f"gpt_result/{config2['try_count']}_50.txt", "w", encoding="utf-8") as f:
    f.write(str(result))

print("查询结果已保存到 query_result.txt 文件中")
