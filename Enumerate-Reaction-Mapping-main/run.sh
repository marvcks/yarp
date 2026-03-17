#!/bin/bash

# 定义基础路径
BASE_INPUT_DIR="yarp-0.5"
BASE_OUTPUT_DIR="all"
SCRIPT_DIR="all"

cd "$SCRIPT_DIR"

# 解析位置参数：
#   ./run.sh N          -> 只运行 lineN
#   ./run.sh START END  -> 运行 lineSTART 到 lineEND
if [ "$#" -eq 1 ]; then
    if ! [[ "$1" =~ ^[0-9]+$ ]]; then
        echo "Error: 参数必须是正整数。"
        echo "Usage: $0 N | $0 START END"
        exit 1
    fi
    START="$1"
    END="$1"
elif [ "$#" -eq 2 ]; then
    if ! [[ "$1" =~ ^[0-9]+$ && "$2" =~ ^[0-9]+$ ]]; then
        echo "Error: 参数必须是正整数。"
        echo "Usage: $0 N | $0 START END"
        exit 1
    fi
    START="$1"
    END="$2"
else
    echo "Usage: $0 N | $0 START END"
    exit 1
fi

if [ "$START" -le 0 ] || [ "$END" -le 0 ]; then
    echo "Error: 参数必须大于 0。"
    exit 1
fi

if [ "$START" -gt "$END" ]; then
    echo "Error: START 不能大于 END。"
    exit 1
fi

for ((i=START; i<=END; i++))
do
    echo "------------------------------------------"
    echo "Processing line${i}..."

    # 定义两种可能的输入文件名
    FILE15="${BASE_INPUT_DIR}/line${i}/reactions_depth15.csv"
    FILE1="${BASE_INPUT_DIR}/line${i}/reactions_depth1.csv"
    OUTPUT_FILE="${BASE_OUTPUT_DIR}/line${i}_aug.csv"

    # 判断哪个文件存在
    if [ -f "$FILE15" ]; then
        INPUT_FILE="$FILE15"
    elif [ -f "$FILE1" ]; then
        INPUT_FILE="$FILE1"
    else
        echo "Warning: No depth1 or depth15 CSV found for line${i}. Skipping."
        continue
    fi

    # 执行 Python 脚本
    echo "Using input: $INPUT_FILE"
    python augment_csv.py \
        "$INPUT_FILE" \
        "$OUTPUT_FILE" \
        --debug \
        --jobs 0

done

echo "All tasks completed."
