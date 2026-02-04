import requests
import json
import os
import uuid


# 定义模型和API配置
api_key = "sk-LifJI8PftxCnfURysfbgcStJVVupWzhnYtKXYBq0x5CcsqDw"

# 定义请求头
headers = {
    "Authorization": f"Bearer {api_key}",
    "Content-Type": "application/json"
}

PROMPT = (
    "In an alkaline environment for the urea oxidation reaction, design a high-entropy "
    "Prussian blue analog catalyst composed of five transition metals. This catalyst should "
    "not only exhibit high catalytic activity and stability but also meet multiple criteria "
    "including low cost, ease of synthesis, environmental friendliness, and corrosion resistance. "
    "Among common metal combinations, considering the balance between catalytic performance and "
    "structural stability, and aiming to avoid certain metals that may pose higher costs or "
    "environmental concerns, a tentative introduction of a metal element with stable valence "
    "states and environmentally benign characteristics is proposed to optimize catalyst performance. "
    "Please analyze in detail the rationale behind the selection of each of the five transition "
    "metal elements, elaborating on their contributions in terms of catalytic activity, stability, "
    "synthetic feasibility, cost-effectiveness, and structural features"
)

API_URL = "https://api.vectorengine.ai/v1/chat/completions"

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

        out_path = os.path.join(out_dir, f"{model}_output_{i}.txt")
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
    ("gpt-5", "gpt")
]

for i in [2,7]:
    for model, out_dir in models:
        call_and_save(model=model, i=i, out_dir=out_dir, headers=headers)

