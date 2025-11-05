#!/usr/bin/env julia

"""
compare_performance.jl

比較階段一和階段二優化的效能差異
"""

using Pkg
Pkg.activate(@__DIR__)

using JSON2MAF
using Printf

function test_stage1(filepath::String, config::FilterConfig)
    println("=" ^ 70)
    println("階段一：預過濾 + 批次寫入")
    println("=" ^ 70)

    GC.gc()
    start_mem = @allocated begin
        start_time = time()

        # 使用階段一優化
        data = parse_nirvana_with_prefilter(filepath, config)

        parse_time = time() - start_time
        variants_count = length(data.positions)

        # 模擬處理
        process_start = time()
        included = 0
        for variant in data.positions
            quality_result = apply_quality_filters(variant, config)
            if !quality_result.passes_quality
                continue
            end

            clinvar = assess_clinvar_pathogenicity(variant.clinvar)
            predictive = assess_predictive_scores(variant, config)
            decision = make_filter_decision(variant, clinvar, predictive)

            if decision.should_include
                included += 1
            end
        end
        process_time = time() - process_start

        total_time = parse_time + process_time

        println(@sprintf("  解析時間:     %.2f 秒", parse_time))
        println(@sprintf("  處理時間:     %.2f 秒", process_time))
        println(@sprintf("  總時間:       %.2f 秒", total_time))
        println(@sprintf("  記憶體分配:   %.2f MB", start_mem / 1024^2))
        println(@sprintf("  解析變異數:   %d", variants_count))
        println(@sprintf("  納入變異數:   %d", included))
        println()

        return (total_time, start_mem, variants_count, included)
    end
end

function test_stage2(filepath::String, config::FilterConfig)
    println("=" ^ 70)
    println("階段二：串流處理")
    println("=" ^ 70)

    GC.gc()
    start_mem = @allocated begin
        start_time = time()

        # 使用階段二優化
        included = 0
        processed = 0

        header = process_nirvana_streaming(filepath, config, batch_size=100) do variant, idx
            processed += 1

            quality_result = apply_quality_filters(variant, config)
            if !quality_result.passes_quality
                return
            end

            clinvar = assess_clinvar_pathogenicity(variant.clinvar)
            predictive = assess_predictive_scores(variant, config)
            decision = make_filter_decision(variant, clinvar, predictive)

            if decision.should_include
                included += 1
            end
        end

        total_time = time() - start_time

        println(@sprintf("  總時間:       %.2f 秒", total_time))
        println(@sprintf("  記憶體分配:   %.2f MB", start_mem / 1024^2))
        println(@sprintf("  處理變異數:   %d", processed))
        println(@sprintf("  納入變異數:   %d", included))
        println()

        return (total_time, start_mem, processed, included)
    end
end

function compare(filepath::String)
    println()
    println("╔" * "═"^68 * "╗")
    println("║" * " "^20 * "效能比較測試" * " "^36 * "║")
    println("╚" * "═"^68 * "╝")
    println()
    println("測試檔案: $filepath")
    println()

    config = create_filter_config()

    # 測試階段一
    (time1, mem1, count1, included1) = test_stage1(filepath, config)

    # 測試階段二
    (time2, mem2, count2, included2) = test_stage2(filepath, config)

    # 比較結果
    println("=" ^ 70)
    println("比較結果")
    println("=" ^ 70)
    println()

    time_improvement = time1 / time2
    mem_improvement = (mem1 - mem2) / mem1 * 100

    println(@sprintf("  速度提升:     %.2fx  (階段二比階段一快 %.1f%%)", time_improvement, (time_improvement - 1) * 100))
    println(@sprintf("  記憶體減少:   %.1f%%", mem_improvement))
    println()

    if included1 == included2
        println("  ✓ 結果一致：兩個版本納入相同數量的變異")
    else
        println("  ⚠ 警告：結果不一致！")
        println(@sprintf("    階段一: %d 個變異", included1))
        println(@sprintf("    階段二: %d 個變異", included2))
    end

    println()
    println("=" ^ 70)
    println("建議")
    println("=" ^ 70)
    println()

    if time_improvement > 1.5
        println("  🎉 階段二優化效果顯著！建議使用串流處理版本。")
    elseif time_improvement > 1.2
        println("  ✓ 階段二有適度改善。對於大型檔案建議使用串流處理。")
    else
        println("  ℹ 階段二改善有限。階段一優化已經足夠。")
    end

    if mem_improvement > 50
        println("  💾 記憶體使用大幅減少！適合處理超大型檔案。")
    elseif mem_improvement > 20
        println("  💾 記憶體使用有明顯改善。")
    end

    println()
end

# 執行比較
if length(ARGS) >= 1
    compare(ARGS[1])
else
    println("用法: julia compare_performance.jl <input.json.gz>")
    println()
    println("範例:")
    println("  julia --project=. compare_performance.jl P.hard-filtered.vcf.annotated.json.gz")
end
