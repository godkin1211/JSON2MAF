"""
parser_tests.jl

Test Nirvana JSON parsing functionality

NOTE: The original synchronous parsing functions (parse_nirvana_json, parse_nirvana_with_prefilter,
process_nirvana_streaming, process_nirvana_parallel, process_nirvana_parallel_no_prefilter,
process_nirvana_streaming_parallel) have been removed from the codebase.

These functions were replaced with process_nirvana_parallel_with_stats() which is tested through
integration tests. The tests below have been disabled as the original functions no longer exist.

If you need to test the new parallel processing functionality, please refer to the integration
tests or examples directory.
"""

using Test
using JSON2MAF

@testset "Nirvana Parser Tests" begin
    # Use relative path from project root
    project_root = dirname(@__DIR__)
    test_file = joinpath(project_root, "examples", "sample_data", "test_mini.json.gz")

    @testset "File Reading Tests" begin
        @test isfile(test_file)

        # NOTE: Tests for parse_nirvana_json() have been removed because this function
        # no longer exists in the codebase. It was replaced with process_nirvana_parallel_with_stats()
        # which uses a different API and is tested through integration tests.

        @test_skip "parse_nirvana_json no longer exists - replaced with process_nirvana_parallel_with_stats()"
    end

    # All remaining test sets have been disabled because they test functions that no longer exist
    # in the codebase (parse_nirvana_json and related synchronous parsing functions).
    #
    # The new parsing architecture uses process_nirvana_parallel_with_stats() which:
    # - Processes files in parallel
    # - Returns statistics instead of full data structures
    # - Writes directly to MAF files
    #
    # For testing the new architecture, see:
    # - Integration tests in other test files
    # - Examples in the examples/ directory

    @testset "Legacy Parser Tests (Disabled)" begin
        @test_skip "Complete parsing test - parse_nirvana_json() no longer exists"
        @test_skip "Header parsing test - parse_nirvana_json() no longer exists"
        @test_skip "Position parsing test - parse_nirvana_json() no longer exists"
        @test_skip "ClinVar parsing test - parse_nirvana_json() no longer exists"
        @test_skip "COSMIC parsing test - parse_nirvana_json() no longer exists"
        @test_skip "Transcript parsing test - parse_nirvana_json() no longer exists"
        @test_skip "Population frequency parsing test - parse_nirvana_json() no longer exists"
        @test_skip "Predictive scores parsing test - parse_nirvana_json() no longer exists"
        @test_skip "dbSNP parsing test - parse_nirvana_json() no longer exists"
        @test_skip "Missing field handling test - parse_nirvana_json() no longer exists"
    end
end

# Real file parsing tests section
@testset "Real File Parsing Tests (Optional)" begin
    project_root = dirname(@__DIR__)
    real_file = joinpath(project_root, "P.hard-filtered.vcf.annotated.json.gz")

    test_real_files = get(ENV, "TEST_REAL_FILES", "0") == "1"

    if isfile(real_file) && test_real_files
        @testset "Real File Basic Parsing" begin
            # NOTE: This test set is disabled because parse_nirvana_json() no longer exists.
            # To process real files, use process_nirvana_parallel_with_stats() instead.
            @test_skip "Parsing large real file - parse_nirvana_json() no longer exists"
        end
    else
        if !isfile(real_file)
            @test_skip "Real file does not exist, skipping test"
        else
            @test_skip "Real file testing disabled (set TEST_REAL_FILES=1 to enable)"
        end
    end
end
