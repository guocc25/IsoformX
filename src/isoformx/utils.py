"""Utility functions for isoform identification."""

import numpy as np
from typing import List, Tuple, Optional, Union
from scipy import stats
import math


def weighted_median(values: List[float], weights: Optional[List[float]] = None) -> float:
    """
    Compute weighted median of values.
    
    Args:
        values: List of values
        weights: Optional list of weights (defaults to uniform weights)
    
    Returns:
        Weighted median value
    """
    if not values:
        raise ValueError("Cannot compute median of empty list")
    
    if weights is None:
        weights = [1.0] * len(values)
    
    if len(values) != len(weights):
        raise ValueError("Values and weights must have same length")
    
    # Sort by values
    sorted_pairs = sorted(zip(values, weights))
    values_sorted, weights_sorted = zip(*sorted_pairs)
    
    # Compute cumulative weights
    cumsum = np.cumsum(weights_sorted)
    total_weight = cumsum[-1]
    
    # Find median position
    median_pos = total_weight / 2
    
    # Find the value at median position
    for i, cum_weight in enumerate(cumsum):
        if cum_weight >= median_pos:
            return values_sorted[i]
    
    return values_sorted[-1]


def jaccard_distance(set1: set, set2: set) -> float:
    """
    Compute Jaccard distance between two sets.
    
    Args:
        set1: First set
        set2: Second set
    
    Returns:
        Jaccard distance (1 - Jaccard similarity)
    """
    if not set1 and not set2:
        return 0.0
    
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    
    if union == 0:
        return 1.0
    
    return 1.0 - (intersection / union)


def edit_distance(seq1: str, seq2: str) -> int:
    """
    Compute edit distance between two sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
    
    Returns:
        Edit distance
    """
    m, n = len(seq1), len(seq2)
    
    # Create distance matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(
                    dp[i - 1][j],      # deletion
                    dp[i][j - 1],      # insertion
                    dp[i - 1][j - 1]   # substitution
                )
    
    return dp[m][n]


def splice_site_distance(sites1: List[int], sites2: List[int], tolerance: int = 10) -> float:
    """
    Compute distance between two sets of splice sites.
    
    Args:
        sites1: First set of splice site positions
        sites2: Second set of splice site positions
        tolerance: Tolerance for matching sites (bp)
    
    Returns:
        Normalized distance (0-1)
    """
    if not sites1 and not sites2:
        return 0.0
    
    if not sites1 or not sites2:
        return 1.0
    
    # Sort sites
    sites1_sorted = sorted(sites1)
    sites2_sorted = sorted(sites2)
    
    # Find matches within tolerance
    matches = 0
    used_indices = set()
    
    for site1 in sites1_sorted:
        for i, site2 in enumerate(sites2_sorted):
            if i in used_indices:
                continue
            if abs(site1 - site2) <= tolerance:
                matches += 1
                used_indices.add(i)
                break
    
    # Compute distance as 1 - (matches / max_sites)
    max_sites = max(len(sites1), len(sites2))
    return 1.0 - (matches / max_sites)


def overlap_coefficient(set1: set, set2: set) -> float:
    """
    Compute overlap coefficient between two sets.
    
    Args:
        set1: First set
        set2: Second set
    
    Returns:
        Overlap coefficient (intersection / min_size)
    """
    if not set1 and not set2:
        return 1.0
    
    intersection = len(set1 & set2)
    min_size = min(len(set1), len(set2))
    
    if min_size == 0:
        return 0.0
    
    return intersection / min_size


def normalize_coverage(coverage: float, read_length: int, transcript_length: int) -> float:
    """
    Normalize coverage by read length and transcript length.
    
    Args:
        coverage: Raw coverage value
        read_length: Length of the read
        transcript_length: Length of the transcript
    
    Returns:
        Normalized coverage
    """
    if transcript_length == 0:
        return 0.0
    
    # Normalize by transcript length and read length
    normalized = coverage * read_length / transcript_length
    return min(normalized, 1.0)  # Cap at 1.0


def compute_confidence_score(
    support_reads: int,
    coverage: float,
    quality_scores: List[float],
    splice_site_consensus: float
) -> float:
    """
    Compute confidence score for an isoform.
    
    Args:
        support_reads: Number of supporting reads
        coverage: Coverage of the isoform
        quality_scores: List of quality scores for supporting reads
        splice_site_consensus: Consensus score for splice sites
    
    Returns:
        Confidence score (0-1)
    """
    if support_reads == 0:
        return 0.0
    
    # Read support component (log scale)
    read_score = min(math.log(support_reads + 1) / math.log(10), 1.0)
    
    # Coverage component
    coverage_score = min(coverage, 1.0)
    
    # Quality component
    avg_quality = np.mean(quality_scores) if quality_scores else 0.0
    quality_score = min(avg_quality, 1.0)
    
    # Splice site consensus component
    consensus_score = splice_site_consensus
    
    # Weighted combination
    weights = [0.3, 0.25, 0.25, 0.2]  # read, coverage, quality, consensus
    scores = [read_score, coverage_score, quality_score, consensus_score]
    
    return sum(w * s for w, s in zip(weights, scores))


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping intervals.
    
    Args:
        intervals: List of (start, end) tuples
    
    Returns:
        List of merged intervals
    """
    if not intervals:
        return []
    
    # Sort by start position
    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]
    
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        
        if start <= last_end:  # Overlapping
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    
    return merged


def compute_exon_overlap(exons1: List[Tuple[int, int]], exons2: List[Tuple[int, int]]) -> float:
    """
    Compute overlap between two sets of exons.
    
    Args:
        exons1: First set of (start, end) tuples
        exons2: Second set of (start, end) tuples
    
    Returns:
        Overlap ratio (0-1)
    """
    if not exons1 or not exons2:
        return 0.0
    
    # Merge intervals for each set
    merged1 = merge_intervals(exons1)
    merged2 = merge_intervals(exons2)
    
    # Compute total overlap
    total_overlap = 0
    total_length1 = sum(end - start + 1 for start, end in merged1)
    total_length2 = sum(end - start + 1 for start, end in merged2)
    
    for start1, end1 in merged1:
        for start2, end2 in merged2:
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            
            if overlap_start <= overlap_end:
                total_overlap += overlap_end - overlap_start + 1
    
    # Return overlap ratio
    max_length = max(total_length1, total_length2)
    if max_length == 0:
        return 0.0
    
    return total_overlap / max_length


def bootstrap_confidence(
    values: List[float],
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95
) -> Tuple[float, float]:
    """
    Compute bootstrap confidence interval for a statistic.
    
    Args:
        values: List of values to bootstrap
        n_bootstrap: Number of bootstrap samples
        confidence_level: Confidence level (e.g., 0.95 for 95%)
    
    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    if not values:
        return (0.0, 0.0)
    
    bootstrap_means = []
    n = len(values)
    
    for _ in range(n_bootstrap):
        # Sample with replacement
        sample = np.random.choice(values, size=n, replace=True)
        bootstrap_means.append(np.mean(sample))
    
    # Compute confidence interval
    alpha = 1 - confidence_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100
    
    lower_bound = np.percentile(bootstrap_means, lower_percentile)
    upper_bound = np.percentile(bootstrap_means, upper_percentile)
    
    return (lower_bound, upper_bound)
