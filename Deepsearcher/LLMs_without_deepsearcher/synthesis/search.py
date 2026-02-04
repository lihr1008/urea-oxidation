import requests
import json
import os
import uuid


# 定义模型和API配置
api_key = "sk-OC6Ph7dHDZdjTkTYgurWptyPYHCKx75PxFr1RxEktnH6Vtxq"

# 定义请求头
headers = {
    "Authorization": f"Bearer {api_key}",
    "Content-Type": "application/json"
}

PROMPT = (
    "How to design a synthesis protocol for a high-entropy Prussian blue analog (HE-PBA) catalyst using cobalt cyanide K3[Co(CN)6] as the framework and incorporating five transition metals (Fe, Co, Ni, Mn, Zn)? Please provide step-by-step instructions."
)

API_URL = "https://api.probex.top/v1/chat/completions"

def call_and_save(model: str, i: int, out_dir: str, headers: dict, prompt: str = PROMPT) -> None:
    os.makedirs(out_dir, exist_ok=True)

    data = {
        "model": model,
        "messages": [{"role": "user", "content": PROMPT}],
        "stream": False,
        "user": str(uuid.uuid4())  # 每次都不同，相当于新会话
    }

    try:
        resp = requests.post(API_URL, headers=headers, json=data, timeout=300)
        if resp.status_code != 200:
            print(f"[{model} #{i}] 请求失败 {resp.status_code}: {resp.text}")
            return

        result = resp.json()
        content = result["choices"][0]["message"]["content"]

        out_path = os.path.join(out_dir, f"output_{i}.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(content)

        print(f"[{model} #{i}] 已写入: {out_path}")

    except requests.exceptions.Timeout:
        print(f"[{model} #{i}] 请求超时：服务器没有响应")
    except requests.exceptions.RequestException as e:
        print(f"[{model} #{i}] 请求发生错误: {e}")
    except (KeyError, ValueError) as e:
        # KeyError: 返回结构不符合预期；ValueError: JSON解析问题
        print(f"[{model} #{i}] 响应解析失败: {e}")

models = [
    ("deepseek-reasoner", "deepseek_result"),
    ("Qwen3-Max", "qwen_result"),
    ("gemini-3-pro-preview", "gemini_result"),
]

for i in [3,4,5]:
    for model, out_dir in models:
        call_and_save(model=model, i=i, out_dir=out_dir, headers=headers)

