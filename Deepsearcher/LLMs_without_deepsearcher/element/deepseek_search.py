import requests
import json

# 定义模型和API配置
model = "deepseek-reasoner"  # 使用的模型，其他模型名请直接复制https://api.probex.top/pricing页面对应的模型名
api_key = "sk-OC6Ph7dHDZdjTkTYgurWptyPYHCKx75PxFr1RxEktnH6Vtxq"

# 定义请求头
headers = {
    "Authorization": f"Bearer {api_key}",
    "Content-Type": "application/json"
}

# 定义请求数据（JSON格式）
data = {
    "model": model,  # 使用前面定义的模型
    "messages": [
        {
            "role": "user",  # 用户发送的消息
            "content": "In an alkaline environment for the urea oxidation reaction, design a high-entropy Prussian blue analog catalyst composed of five transition metals. This catalyst should not only exhibit high catalytic activity and stability but also meet multiple criteria including low cost, ease of synthesis, environmental friendliness, and corrosion resistance. Among common metal combinations, considering the balance between catalytic performance and structural stability, and aiming to avoid certain metals that may pose higher costs or environmental concerns, a tentative introduction of a metal element with stable valence states and environmentally benign characteristics is proposed to optimize catalyst performance. Please analyze in detail the rationale behind the selection of each of the five transition metal elements, elaborating on their contributions in terms of catalytic activity, stability, synthetic feasibility, cost-effectiveness, and structural features"  # 消息内容
        }
    ],
    "stream": False  # 是否启用流式响应（False 表示一次性返回完整数据，True 表示流式返回数据）
}

# 发送POST请求
try:
    response = requests.post(
        "https://api.probex.top/v1/chat/completions",
        headers=headers,
        json=data,
        timeout=300
    )

    # 处理非流式响应
    if response.status_code == 200:
        result = response.json()  # 解析JSON响应
        content = result['choices'][0]['message']['content']  # 提取content字段
        print(content)  # 直接打印content内容
        out_path = "deepseek/deepseek_output_11.txt"
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"已写入: {out_path}")

    else:
        print("请求失败，状态码:", response.status_code)
        print("错误信息:", response.text)

except requests.exceptions.Timeout:
    print("请求超时：服务器没有响应")
except requests.exceptions.RequestException as e:
    print(f"请求发生错误: {e}")


