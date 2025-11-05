"""
MAFWriter.jl

MAF 檔案寫入模組
將 MAFRecord 寫入為標準 MAF 格式檔案 (TSV)
"""

module MAFWriter

using CSV
using DataFrames
using ..DataStructures

export write_maf_file, maf_header, mafrecord_to_row, create_maf_writer, write_maf_batch, close_maf_writer, get_total_written, merge_maf_files, count_maf_records

"""
MAF 檔案標準欄位順序
"""
const MAF_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "COSMIC_ID",
    "ClinVar_ID",
    "ClinVar_Review_Status",
    "ClinVar_Significance",
    "ClinVar_Disease",
    "PrimateAI_Score",
    "DANN_Score",
    "REVEL_Score",
    "gnomAD_AF",
    "gnomAD_EAS_AF",
    "Depth",
    "VAF"
]

"""
    maf_header() -> Vector{String}

返回 MAF 檔案的標準表頭
"""
function maf_header()
    return MAF_COLUMNS
end

"""
    format_maf_value(val) -> String

格式化 MAF 欄位值
- nothing 或空字串 → "."
- 其他值 → 轉為字串
"""
function format_maf_value(val)
    if val === nothing || val == ""
        return "."
    else
        return string(val)
    end
end

"""
    mafrecord_to_row(record::MAFRecord) -> Vector{String}

將 MAFRecord 轉換為 MAF 檔案的一行資料

# 參數
- `record`: MAFRecord 物件

# 返回
- `Vector{String}`: 按照 MAF_COLUMNS 順序排列的欄位值
"""
function mafrecord_to_row(record::MAFRecord)::Vector{String}
    return [
        format_maf_value(record.hugo_symbol),
        format_maf_value(record.chromosome),
        format_maf_value(record.start_position),
        format_maf_value(record.end_position),
        format_maf_value(record.strand),
        format_maf_value(record.variant_classification),
        format_maf_value(record.variant_type),
        format_maf_value(record.reference_allele),
        format_maf_value(record.tumor_seq_allele1),
        format_maf_value(record.tumor_seq_allele2),
        format_maf_value(record.tumor_sample_barcode),
        format_maf_value(record.hgvsc),
        format_maf_value(record.hgvsp),
        format_maf_value(record.hgvsp_short),
        format_maf_value(record.transcript_id),
        format_maf_value(record.dbsnp_rs),
        format_maf_value(record.dbsnp_val_status),
        format_maf_value(record.cosmic_id),
        format_maf_value(record.clinvar_id),
        format_maf_value(record.clinvar_review_status),
        format_maf_value(record.clinvar_significance),
        format_maf_value(record.clinvar_disease),
        format_maf_value(record.primate_ai_score),
        format_maf_value(record.dann_score),
        format_maf_value(record.revel_score),
        format_maf_value(record.gnomad_af),
        format_maf_value(record.gnomad_eas_af),
        format_maf_value(record.depth),
        format_maf_value(record.vaf)
    ]
end

"""
    write_maf_file(output_path::String, records::Vector{MAFRecord}; append::Bool=false)

將 MAFRecord 列表寫入 MAF 檔案

# 參數
- `output_path`: 輸出檔案路徑
- `records`: MAFRecord 向量
- `append`: 是否追加模式（預設 false，會覆蓋現有檔案）

# 範例
```julia
records = [maf_record1, maf_record2, ...]
write_maf_file("output.maf", records)
```

# 注意
- MAF 格式為 TSV (tab-separated values)
- 第一行為標準表頭
- 空值以 "." 表示
"""
function write_maf_file(output_path::String, records::Vector{MAFRecord}; append::Bool=false)
    # 如果沒有記錄，只寫入表頭
    if isempty(records)
        open(output_path, "w") do io
            println(io, join(maf_header(), "\t"))
        end
        return
    end

    # 轉換所有記錄為行資料
    rows = [mafrecord_to_row(record) for record in records]

    # 建立 DataFrame
    df = DataFrame([col => [row[i] for row in rows] for (i, col) in enumerate(MAF_COLUMNS)])

    # 寫入檔案
    if append && isfile(output_path)
        # 追加模式：不寫入表頭
        CSV.write(output_path, df, delim='\t', append=true)
    else
        # 覆蓋模式：寫入表頭
        CSV.write(output_path, df, delim='\t')
    end
end

"""
    write_maf_file(output_path::String, record::MAFRecord; append::Bool=false)

寫入單一 MAFRecord（便利函數）
"""
function write_maf_file(output_path::String, record::MAFRecord; append::Bool=false)
    write_maf_file(output_path, [record], append=append)
end

# ============================================================================
# 批次寫入 API - 用於優化大檔案處理效能
# ============================================================================

"""
    MAFBatchWriter

批次 MAF 寫入器，用於優化大量記錄的寫入效能

# 欄位
- `output_path`: 輸出檔案路徑
- `io`: 檔案句柄
- `buffer`: 記錄緩衝區
- `batch_size`: 批次大小
- `total_written`: 已寫入記錄總數
"""
mutable struct MAFBatchWriter
    output_path::String
    io::Union{IO, Nothing}
    buffer::Vector{MAFRecord}
    batch_size::Int
    total_written::Int
end

"""
    create_maf_writer(output_path::String; batch_size::Int=1000) -> MAFBatchWriter

創建批次 MAF 寫入器

# 參數
- `output_path`: 輸出檔案路徑
- `batch_size`: 批次大小（預設 1000），累積這麼多記錄後才寫入

# 返回
- `MAFBatchWriter`: 批次寫入器實例

# 使用方式
```julia
writer = create_maf_writer("output.maf", batch_size=1000)
for record in records
    write_maf_batch(writer, record)
end
close_maf_writer(writer)  # 別忘了關閉！
```

# 效能優勢
- 減少 I/O 次數
- 批次處理提升效率
- 適合處理大量記錄
"""
function create_maf_writer(output_path::String; batch_size::Int=1000)::MAFBatchWriter
    # 打開檔案並寫入表頭
    io = open(output_path, "w")
    println(io, join(maf_header(), "\t"))

    return MAFBatchWriter(
        output_path,
        io,
        MAFRecord[],
        batch_size,
        0
    )
end

"""
    write_maf_batch(writer::MAFBatchWriter, record::MAFRecord)

將記錄加入批次寫入器

當緩衝區達到 batch_size 時自動寫入檔案。

# 參數
- `writer`: 批次寫入器
- `record`: 要寫入的 MAF 記錄
"""
function write_maf_batch(writer::MAFBatchWriter, record::MAFRecord)
    push!(writer.buffer, record)

    # 當緩衝區滿時，寫入檔案
    if length(writer.buffer) >= writer.batch_size
        flush_buffer(writer)
    end
end

"""
    flush_buffer(writer::MAFBatchWriter)

將緩衝區的記錄寫入檔案
"""
function flush_buffer(writer::MAFBatchWriter)
    if isempty(writer.buffer)
        return
    end

    # 將所有緩衝的記錄轉換為行並寫入
    for record in writer.buffer
        row = mafrecord_to_row(record)
        println(writer.io, join(row, "\t"))
    end

    writer.total_written += length(writer.buffer)
    empty!(writer.buffer)  # 清空緩衝區
end

"""
    close_maf_writer(writer::MAFBatchWriter)

關閉批次寫入器並確保所有資料已寫入

# 參數
- `writer`: 批次寫入器

# 注意
必須呼叫此函數以確保所有緩衝的記錄都寫入檔案！
"""
function close_maf_writer(writer::MAFBatchWriter)
    # 寫入剩餘的緩衝記錄
    flush_buffer(writer)

    # 關閉檔案
    if writer.io !== nothing
        close(writer.io)
        writer.io = nothing
    end
end

"""
    get_total_written(writer::MAFBatchWriter) -> Int

返回已寫入的記錄總數
"""
function get_total_written(writer::MAFBatchWriter)::Int
    return writer.total_written + length(writer.buffer)
end

# ============================================================================
# MAF 檔案合併 - 用於多執行緒並行處理
# ============================================================================

"""
    merge_maf_files(input_files::Vector{String}, output_file::String; keep_temp::Bool=false)

合併多個 MAF 檔案為單一檔案

這個函數用於合併多執行緒並行處理產生的多個 MAF 檔案。
會保留所有記錄，並確保只有一個標準表頭。

# 參數
- `input_files`: 輸入 MAF 檔案路徑列表
- `output_file`: 輸出合併後的 MAF 檔案路徑
- `keep_temp`: 是否保留臨時檔案（預設 false，會刪除）

# 使用範例
```julia
# 並行處理後合併
merge_maf_files(
    ["output_1.maf", "output_2.maf", "output_3.maf"],
    "final_output.maf",
    keep_temp=false
)
```

# 注意事項
- 輸入檔案必須都存在
- 所有輸入檔案必須有相同的欄位結構
- 預設會刪除臨時檔案以節省空間
"""
function merge_maf_files(input_files::Vector{String}, output_file::String; keep_temp::Bool=false)
    # 檢查輸入檔案
    for file in input_files
        if !isfile(file)
            error("Input file not found: $file")
        end
    end

    # 開啟輸出檔案並寫入表頭
    open(output_file, "w") do out_io
        # 寫入標準表頭
        println(out_io, join(maf_header(), "\t"))

        # 逐一讀取並合併每個輸入檔案
        for (idx, input_file) in enumerate(input_files)
            open(input_file, "r") do in_io
                # 跳過輸入檔案的表頭
                readline(in_io)

                # 複製所有資料行
                while !eof(in_io)
                    line = readline(in_io)
                    if !isempty(strip(line))
                        println(out_io, line)
                    end
                end
            end
        end
    end

    # 刪除臨時檔案（如果需要）
    if !keep_temp
        for file in input_files
            try
                rm(file)
            catch e
                @warn "Failed to delete temporary file: $file" exception=e
            end
        end
    end

    return output_file
end

"""
    count_maf_records(filepath::String) -> Int

計算 MAF 檔案中的記錄數（不含表頭）

# 參數
- `filepath`: MAF 檔案路徑

# 返回
- `Int`: 記錄數量
"""
function count_maf_records(filepath::String)::Int
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    count = 0
    open(filepath, "r") do io
        # 跳過表頭
        readline(io)

        # 計數非空行
        while !eof(io)
            line = readline(io)
            if !isempty(strip(line))
                count += 1
            end
        end
    end

    return count
end

end # module MAFWriter
