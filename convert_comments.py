#!/usr/bin/env python3
"""
Script to convert Chinese comments to English in Julia files
"""

import re
import sys

# Translation mappings for common terms
translations = {
    # General terms
    "模組": "module",
    "解析": "parse/parsing",
    "檔案": "file",
    "支援": "supports",
    "高效": "efficient",
    "參數": "Parameters",
    "返回": "Returns",
    "範例": "Example",
    "效能優勢": "Performance Benefits",
    "注意事項": "Notes",
    "使用範例": "Usage Example",
    "使用方式": "Usage",
    "欄位": "Fields/field",
    "函數": "function",
    "變異": "variant",
    "位點": "position",
    "資料": "data",
    "結構": "structure",
    "配置": "configuration",
    "過濾": "filter/filtering",
    "評估": "assess/assessment",
    "致病性": "pathogenicity",
    "分數": "score",
    "閾值": "threshold",
    "品質": "quality",
    "頻率": "frequency",
    "族群": "population",
    "東亞": "East Asian",
    "測序": "sequencing",
    "深度": "depth",
    "寫入": "write/writing",
    "讀取": "read/reading",
    "建立": "create/build",
    "檢查": "check",
    "判斷": "determine/judge",
    "提取": "extract",
    "格式": "format",
    "轉換": "convert/conversion",
    "合併": "merge",
    "優化": "optimization/optimize",
    "記憶體": "memory",
    "處理": "process/processing",
    "速度": "speed",
    "批次": "batch",
    "並行": "parallel",
    "執行緒": "thread",
    "串流": "streaming",

    # Specific terms
    "Nirvana JSON 檔案解析模組": "Nirvana JSON file parsing module",
    "支援 gzipped JSON 檔案的高效解析": "Supports efficient parsing of gzipped JSON files",
    "解析 Nirvana gzipped JSON 檔案": "Parse Nirvana gzipped JSON file",
    "Nirvana JSON.gz 檔案的路徑": "Path to Nirvana JSON.gz file",
    "包含 header, positions, genes 的完整資料結構": "Complete data structure containing header, positions, genes",
    "解析 Nirvana JSON 的 header 區塊": "Parse Nirvana JSON header section",
    "解析所有變異位點": "Parse all variant positions",
    "解析單一變異位點": "Parse single variant position",
    "基本資訊": "Basic information",
    "樣本資訊": "Sample information",
    "變異註解": "Variant annotations",
    "轉錄本": "transcript",
    "族群頻率": "population frequency",
    "預測分數": "predictive score",
    "測序品質": "sequencing quality",
    "通過所有初步檢查": "Passed all initial checks",
    "讀取並解析 JSON": "Read and parse JSON",
    "純文字 JSON": "Plain text JSON",
    "解析各區塊": "Parse sections",
    "區塊是可選的": "section is optional",
    "檢查項目": "Check items",
    "可以跳過這個變異": "Can skip this variant",
    "應該進行完整解析": "Should perform full parsing",
    "深度不足，跳過": "Insufficient depth, skip",
    "過低，跳過": "Too low, skip",
    "過高，跳過": "Too high, skip",
    "取第一個樣本的資料": "Get data from first sample",
    "處理第一個變異": "Process first variant",
    "通常只有一個": "Usually only one",
    "沒有變異註解的情況": "Case with no variant annotations",
    "預設": "Default",
    "可能是": "May be",
    "如果包含": "If contains",
    "取第一個數字": "Take first number",
    "如果是範圍": "If is range",
    "直接嘗試解析": "Try parsing directly",
    "嘗試解析為整數": "Try parsing as integer",
    "失敗則返回": "Return on failure",
    "第一個替代等位基因": "First alternate allele",
    "避免與 Julia 關鍵字衝突": "Avoid conflict with Julia keyword",
}

# Specific docstring translations
docstring_translations = {
    "ClinVar 致病性評估與優先級判斷模組": "ClinVar pathogenicity assessment and priority determination module",
    "根據 review status 與 clinical significance 判斷變異的致病性": "Determine variant pathogenicity based on review status and clinical significance",
    "整合過濾決策引擎": "Integrated filtering decision engine",
    "整合 ClinVar 與預測分數評估，做出最終納入/排除決策": "Integrate ClinVar and predictive score assessment to make final include/exclude decisions",
    "預測分數評估模組": "Predictive scores assessment module",
    "整合多種預測分數": "Integrate multiple predictive scores",
    "提供 likely pathogenic 補充判斷": "Provide supplementary likely pathogenic determination",
    "品質與族群頻率預過濾模組": "Quality and population frequency pre-filtering module",
    "在進行 ClinVar 評估前先過濾低品質與常見變異": "Filter low quality and common variants before ClinVar assessment",
    "定義 JSON2MAF 專案的核心資料結構": "Define core data structures for JSON2MAF project",
    "過濾配置管理模組": "Filter configuration management module",
    "提供配置的創建、驗證和顯示功能": "Provides configuration creation, validation and display functions",
    "MAF 檔案寫入模組": "MAF file writing module",
    "將 MAFRecord 寫入為標準 MAF 格式檔案": "Write MAFRecord to standard MAF format file",
}

def main():
    print("Translation mappings loaded")
    print(f"Total basic translations: {len(translations)}")
    print(f"Total docstring translations: {len(docstring_translations)}")
    print("\nReady to translate files manually using the mappings above")

if __name__ == "__main__":
    main()
