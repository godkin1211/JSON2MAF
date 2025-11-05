#!/usr/bin/env julia

"""
效能測試腳本
比較優化前後的處理速度和記憶體使用
"""

using Pkg
Pkg.activate(@__DIR__)

using JSON2MAF
using Printf

function test_performance(filepath::String, output_path::String)
    println("=" ^ 70)
    println("JSON2MAF 效能測試")
    println("=" ^ 70)
    println()

    # 建立配置
    config = create_filter_config()

    println("測試檔案: $filepath")
    println("輸出檔案: $output_path")
    println()

    # 測試優化版本（帶預過濾）
    println("🚀 測試優化版本（預過濾 + 批次寫入）")
    println("-" ^ 70)

    GC.gc()  # 清理垃圾收集
    start_mem = @allocated begin
        start_time = time()

        data = parse_nirvana_with_prefilter(filepath, config)

        parse_time = time() - start_time
        println(@sprintf("  ⏱  解析時間: %.2f 秒", parse_time))
        println(@sprintf("  📊 解析後變異數: %d", length(data.positions)))

        # 計算過濾率
        println()
        println("  開始處理過濾...")

        process_start = time()
        included_count = 0

        for variant in data.positions
            quality_result = apply_quality_filters(variant, config)
            if !quality_result.passes_quality
                continue
            end

            clinvar_assessment = assess_clinvar_pathogenicity(variant.clinvar)
            predictive_assessment = assess_predictive_scores(variant, config)
            decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

            if decision.should_include
                included_count += 1
            end
        end

        process_time = time() - process_start
        println(@sprintf("  ⏱  過濾時間: %.2f 秒", process_time))
        println(@sprintf("  ✅ 最終納入: %d 個變異", included_count))
        println()
        println(@sprintf("  📈 總處理時間: %.2f 秒", parse_time + process_time))
    end

    println(@sprintf("  💾 記憶體分配: %.2f MB", start_mem / 1024^2))
    println()

    println("=" ^ 70)
    println("測試完成！")
    println("=" ^ 70)
end

# 執行測試
if length(ARGS) >= 1
    input_file = ARGS[1]
    output_file = length(ARGS) >= 2 ? ARGS[2] : "test_output.maf"
    test_performance(input_file, output_file)
else
    println("用法: julia benchmark_test.jl <input.json.gz> [output.maf]")
    println()
    println("範例:")
    println("  julia --project=. benchmark_test.jl P.hard-filtered.vcf.annotated.json.gz test_output.maf")
end
