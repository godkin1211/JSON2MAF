#!/usr/bin/env julia

"""
json2maf.jl

Nirvana JSON 致癌性重要位點過濾工具

使用多執行緒並行處理大規模 Nirvana JSON 檔案：
- 利用多核心 CPU 並行處理變異
- 完整的統計追蹤（包括預過濾階段）
- 每個執行緒獨立處理一批變異並合併結果
- 支援大型 VCF 註解檔案的高效能處理

⚠️ 重要：建議設定執行緒數以獲得最佳效能
export JULIA_NUM_THREADS=8

用法：
    julia -t 8 bin/json2maf.jl -i input.json.gz -o output.maf [options]
"""

using Pkg
# 啟用專案環境
Pkg.activate(joinpath(@__DIR__, ".."))

using ArgParse
using JSON2MAF

"""
    parse_commandline()

解析命令行參數
"""
function parse_commandline()
    s = ArgParseSettings(
        prog = "json2maf",
        description = "Nirvana JSON 致癌性重要位點過濾工具",
        version = string(JSON2MAF.version()),
        add_version = true
    )

    @add_arg_table! s begin
        "--input", "-i"
            help = "輸入 Nirvana JSON.gz 檔案路徑"
            required = true
            arg_type = String

        "--output", "-o"
            help = "輸出 MAF 檔案路徑"
            required = true
            arg_type = String

        "--min-depth"
            help = "最小測序深度 (預設: 30)"
            arg_type = Int
            default = 30

        "--min-vaf"
            help = "最小變異等位基因頻率 VAF (預設: 0.03)"
            arg_type = Float64
            default = 0.03

        "--max-eas-af"
            help = "東亞族群最大等位基因頻率 (預設: 0.01)"
            arg_type = Float64
            default = 0.01

        "--min-revel"
            help = "REVEL 分數閾值 (預設: 0.75)"
            arg_type = Float64
            default = 0.75

        "--min-primate-ai"
            help = "PrimateAI-3D 分數閾值 (預設: 0.8)"
            arg_type = Float64
            default = 0.8

        "--min-dann"
            help = "DANN 分數閾值 (預設: 0.96)"
            arg_type = Float64
            default = 0.96

        "--keep-temp"
            help = "保留臨時檔案（預設會刪除）"
            action = :store_true

        "--stats"
            help = "統計報告輸出路徑 (可選)"
            arg_type = String
            default = nothing

        "--verbose", "-v"
            help = "詳細輸出模式"
            action = :store_true

        "--quiet", "-q"
            help = "安靜模式（不顯示進度）"
            action = :store_true
    end

    return parse_args(s)
end

"""
    print_statistics(stats::Dict, output_path::Union{String, Nothing}=nothing)

輸出統計資訊
"""
function print_statistics(stats::Dict, output_path::Union{String, Nothing}=nothing)
    report = """

    ═══════════════════════════════════════════════════════════
                        過濾統計報告
    ═══════════════════════════════════════════════════════════

    處理模式:               多執行緒並行處理
    執行緒數:               $(stats[:num_threads])

    品質過濾:
      - 通過品質過濾:       $(stats[:passed_quality])
      - 深度不足:           $(stats[:failed_depth])
      - VAF 過低:           $(stats[:failed_vaf])
      - 族群頻率過高:       $(stats[:failed_af])

    致病性評估:
      - ClinVar Pathogenic:         $(stats[:clinvar_pathogenic])
      - ClinVar Likely pathogenic:  $(stats[:clinvar_likely])
      - 預測分數支持:               $(stats[:predictive_likely])
        * PrimateAI-3D 單獨支持:    $(stats[:primate_ai_only])
        * 2+ 分數支持:              $(stats[:multi_score])

    最終結果:
      - 納入變異數:         $(stats[:included])
      - 排除變異數:         $(stats[:excluded])

    ═══════════════════════════════════════════════════════════
    """

    println(report)

    # 如果指定了統計檔案路徑，寫入檔案
    if output_path !== nothing
        open(output_path, "w") do io
            write(io, report)
        end
        println("統計報告已寫入: $output_path")
    end
end

"""
    process_nirvana_json(input_path::String, output_path::String, config::FilterConfig;
                        keep_temp::Bool=false, verbose::Bool=false, quiet::Bool=false,
                        stats_path::Union{String,Nothing}=nothing)

處理 Nirvana JSON 檔案並輸出為 MAF 格式
"""
function process_nirvana_json(input_path::String, output_path::String, config::FilterConfig;
                                      keep_temp::Bool=false, verbose::Bool=false, quiet::Bool=false,
                                      stats_path::Union{String,Nothing}=nothing)

    nthreads = Threads.nthreads()

    if verbose
        println("\n開始處理: $input_path")
        println("輸出檔案: $output_path")
        println("處理模式: 多執行緒並行處理")
        println("執行緒數: $nthreads")
        display_config(config)
    end

    if nthreads == 1
        @warn "只有 1 個執行緒！建議設定 JULIA_NUM_THREADS 以獲得更好的效能。"
        @warn "例如: export JULIA_NUM_THREADS=8"
    end

    # 為每個執行緒創建獨立的統計資料
    thread_stats = [Dict(
        :passed_quality => 0,
        :failed_depth => 0,
        :failed_vaf => 0,
        :failed_af => 0,
        :clinvar_pathogenic => 0,
        :clinvar_likely => 0,
        :predictive_likely => 0,
        :primate_ai_only => 0,
        :multi_score => 0,
        :included => 0,
        :excluded => 0
    ) for _ in 1:nthreads]

    # 為每個執行緒創建獨立的 MAF 寫入器
    temp_files = ["$(output_path).thread_$(i).tmp" for i in 1:nthreads]
    thread_writers = [create_maf_writer(temp_files[i], batch_size=1000) for i in 1:nthreads]

    if verbose
        println("\n[1/4] 並行處理 Nirvana JSON...")
    end

    # 並行處理所有變異（在預過濾階段就追蹤完整統計）
    try
        header = process_nirvana_parallel_with_stats(input_path, config, thread_stats) do variant, idx, tid
            stats = thread_stats[tid]
            writer = thread_writers[tid]

            # 品質過濾
            quality_result = apply_quality_filters(variant, config)

            if !quality_result.passes_quality
                # 更新統計 - 根據失敗原因分類
                if quality_result.failure_reason !== nothing
                    reason_lower = lowercase(quality_result.failure_reason)
                    if occursin("depth", reason_lower)
                        stats[:failed_depth] += 1
                    elseif occursin("frequency", reason_lower) || occursin("vaf", reason_lower)
                        stats[:failed_vaf] += 1
                    elseif occursin("population", reason_lower) ||
                           occursin("allele frequency", reason_lower) ||
                           occursin("east asian", reason_lower)
                        stats[:failed_af] += 1
                    end
                end
                return
            end

            stats[:passed_quality] += 1

            # ClinVar 評估
            clinvar_assessment = assess_clinvar_pathogenicity(variant.clinvar)

            # 預測分數評估
            predictive_assessment = assess_predictive_scores(variant, config)

            # 整合決策
            decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

            # 更新統計並寫入
            if decision.should_include
                stats[:included] += 1

                if decision.primary_evidence == "ClinVar"
                    if decision.pathogenicity_class == "Pathogenic"
                        stats[:clinvar_pathogenic] += 1
                    elseif decision.pathogenicity_class == "Likely pathogenic"
                        stats[:clinvar_likely] += 1
                    end
                elseif decision.primary_evidence == "Predictive"
                    stats[:predictive_likely] += 1
                    if has_primate_ai_support(predictive_assessment) &&
                       count_supporting_predictive_scores(predictive_assessment) == 1
                        stats[:primate_ai_only] += 1
                    else
                        stats[:multi_score] += 1
                    end
                end

                # 轉換為 MAF 記錄並寫入
                maf_record = variant_to_maf(variant, decision)
                write_maf_batch(writer, maf_record)

                # 批次輸出（每 50 個變異輸出一次，減少 I/O 開銷）
                if verbose && stats[:included] % 50 == 0
                    println("  ✓ [Thread $tid] Processed $(stats[:included]) variants " *
                           "(latest: $(maf_record.hugo_symbol) $(maf_record.chromosome):$(maf_record.start_position))")
                end
            else
                stats[:excluded] += 1
            end
        end

    finally
        # 關閉所有寫入器
        for writer in thread_writers
            close_maf_writer(writer)
        end
    end

    if verbose
        println("\n[2/4] 合併執行緒輸出...")
    end

    # 合併所有執行緒的輸出檔案
    merge_maf_files(temp_files, output_path, keep_temp=keep_temp)

    # 合併統計資料
    total_stats = Dict(
        :num_threads => nthreads,
        :passed_quality => sum(s[:passed_quality] for s in thread_stats),
        :failed_depth => sum(s[:failed_depth] for s in thread_stats),
        :failed_vaf => sum(s[:failed_vaf] for s in thread_stats),
        :failed_af => sum(s[:failed_af] for s in thread_stats),
        :clinvar_pathogenic => sum(s[:clinvar_pathogenic] for s in thread_stats),
        :clinvar_likely => sum(s[:clinvar_likely] for s in thread_stats),
        :predictive_likely => sum(s[:predictive_likely] for s in thread_stats),
        :primate_ai_only => sum(s[:primate_ai_only] for s in thread_stats),
        :multi_score => sum(s[:multi_score] for s in thread_stats),
        :included => sum(s[:included] for s in thread_stats),
        :excluded => sum(s[:excluded] for s in thread_stats)
    )

    total_written = count_maf_records(output_path)

    if verbose
        println("\n[3/4] 完成 MAF 檔案寫入...")
        println("  ✓ 成功合併並寫入 $total_written 筆記錄到 $output_path")
        if !keep_temp
            println("  ✓ 已刪除臨時檔案")
        end
    end

    # 輸出統計
    if verbose
        println("\n[4/4] 生成統計報告...")
    end

    print_statistics(total_stats, stats_path)

    return total_stats
end

"""
    main()

主程式入口
"""
function main()
    # 解析參數
    args = parse_commandline()

    # 建立配置
    config = create_filter_config(
        min_total_depth = args["min-depth"],
        min_variant_frequency = args["min-vaf"],
        max_eas_af = args["max-eas-af"],
        min_revel_score = args["min-revel"],
        min_primate_ai_score = args["min-primate-ai"],
        min_dann_score = args["min-dann"]
    )

    # 驗證配置
    try
        validate_config(config)
    catch e
        println(stderr, "錯誤: 配置參數無效 - $e")
        exit(1)
    end

    # 檢查輸入檔案
    if !isfile(args["input"])
        println(stderr, "錯誤: 輸入檔案不存在: $(args["input"])")
        exit(1)
    end

    # 檢查執行緒數
    nthreads = Threads.nthreads()
    if nthreads == 1 && !args["quiet"]
        println("警告: 只有 1 個執行緒！")
        println("      建議設定 JULIA_NUM_THREADS 環境變數以啟用多執行緒")
        println("      例如: export JULIA_NUM_THREADS=8")
        println("      或使用: julia -t 8 $(PROGRAM_FILE)")
        println()
    end

    try
        # 處理檔案
        process_nirvana_json(
            args["input"],
            args["output"],
            config,
            keep_temp = args["keep-temp"],
            verbose = args["verbose"],
            quiet = args["quiet"],
            stats_path = args["stats"]
        )

        println("\n✓ 處理完成！（使用 $nthreads 個執行緒並行處理）")
        exit(0)

    catch e
        println(stderr, "\n✗ 發生錯誤: $e")
        if args["verbose"]
            showerror(stderr, e, catch_backtrace())
        end
        exit(1)
    end
end

# 執行主程式
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
