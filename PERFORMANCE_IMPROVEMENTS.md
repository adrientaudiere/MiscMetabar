# Performance Improvements in MiscMetabar

This document summarizes the speed optimizations implemented in MiscMetabar to enhance performance without breaking existing functionality.

## Summary of Optimizations

### 1. Pre-allocated List Initialization (Memory Efficiency)

**Problem**: Using `list()` and growing lists iteratively is inefficient as it requires repeated memory reallocation.

**Solution**: Replace `list()` with pre-allocated `vector("list", length)` when the final size is known.

**Files Modified**:
- `R/alpha_div_test.R`: Pre-allocated `res_perm` and `p_perm` lists
- `R/dada_phyloseq.R`: Pre-allocated `cmd` list based on file count
- `R/plot_functions.R`: Pre-allocated `p` list for biplot pairs

**Performance Gain**: Reduces memory fragmentation and eliminates repeated reallocation overhead.

### 2. Efficient Data Structure Building (Vectorization)

**Problem**: Repeated `rbind()` operations in loops have O(n²) complexity due to copying.

**Solution**: Collect results in a list, then use single `do.call(rbind, list)` operation.

**Files Modified**:
- `R/plot_functions.R`: Optimized sankey diagram data preparation loops

**Performance Gain**: Changes from O(n²) to O(n) complexity for data frame building.

### 3. Vectorized Matrix Operations

**Problem**: Nested loops for matrix operations are slow and inefficient.

**Solution**: Replace nested loops with vectorized operations using matrix indexing.

**Files Modified**:
- `R/plot_functions.R`: Completely vectorized `net_matrix2links` function

**Performance Gain**: Eliminates nested loops, uses efficient vectorized operations.

### 4. Optimized Loop Patterns

**Problem**: Using `seq_len(length(x))` instead of `seq_along(x)` and `seq(1, length(x))`.

**Solution**: Replace with more efficient and readable `seq_along(x)` pattern.

**Files Modified**:
- `R/beta_div_test.R`: Fixed loop iterators for sam_variables
- `R/plot_functions.R`: Fixed multiple loop patterns
- `R/dada_phyloseq.R`: Fixed guild processing loop

**Performance Gain**: Slightly more efficient and more readable code.

### 5. Reduced Redundant Computations

**Problem**: Computing the same `sapply()` result multiple times.

**Solution**: Cache computation results and reuse them.

**Files Modified**:
- `R/beta_div_test.R`: Pre-computed sapply results for R-square calculations

**Performance Gain**: Eliminates redundant expensive computations.

### 6. Optimized String Operations

**Problem**: Nested `gsub()` calls for simple pattern replacements.

**Solution**: Use single regex pattern to replace multiple variations efficiently.

**Files Modified**:
- `R/dada_phyloseq.R`: Replaced nested gsub with single `sub()` pattern for file extensions

**Performance Gain**: Single regex operation instead of multiple string passes.

## Testing

All optimizations were tested with `test_performance_improvements.R` to ensure:
- Functionality is preserved
- Results are identical to original implementations
- Edge cases are handled correctly

## Impact

These optimizations provide:

1. **Memory Efficiency**: Reduced memory allocation overhead
2. **Computational Efficiency**: Faster loops and vectorized operations  
3. **Scalability**: Better performance with larger datasets
4. **Maintainability**: More readable and standard R idioms

## Backward Compatibility

All changes maintain 100% backward compatibility:
- Function signatures unchanged
- Return values identical
- No breaking changes to user-facing APIs

The optimizations are purely internal implementation improvements that enhance performance while preserving all existing functionality.