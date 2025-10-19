"""Confidence scoring for isoforms."""

from typing import List, Dict, Optional, Tuple
import numpy as np
from dataclasses import dataclass
from .models import Isoform, ReadSignature
from .utils import compute_confidence_score, bootstrap_confidence


@dataclass
class ScoringParams:
    """Parameters for confidence scoring."""
    min_reads: int = 2
    min_coverage: float = 0.3
    min_quality: float = 0.2
    splice_consensus_weight: float = 0.3
    read_support_weight: float = 0.25
    coverage_weight: float = 0.25
    quality_weight: float = 0.2


class ConfidenceScorer:
    """Computes confidence scores for isoforms."""
    
    def __init__(self, params: ScoringParams):
        self.params = params
    
    def score_isoform(self, isoform: Isoform) -> float:
        """
        Compute confidence score for an isoform.
        
        Args:
            isoform: Isoform to score
        
        Returns:
            Confidence score (0-1)
        """
        if not isoform.read_signatures:
            return 0.0
        
        # Extract metrics
        read_support = len(isoform.read_signatures)
        coverage = isoform.coverage
        quality_scores = [r.quality_score for r in isoform.read_signatures]
        
        # Compute splice site consensus
        splice_consensus = self._compute_splice_consensus(isoform)
        
        # Compute individual component scores
        read_score = self._compute_read_support_score(read_support)
        coverage_score = self._compute_coverage_score(coverage)
        quality_score = self._compute_quality_score(quality_scores)
        consensus_score = self._compute_consensus_score(splice_consensus)
        
        # Weighted combination
        weights = [
            self.params.read_support_weight,
            self.params.coverage_weight,
            self.params.quality_weight,
            self.params.splice_consensus_weight
        ]
        
        scores = [read_score, coverage_score, quality_score, consensus_score]
        
        confidence = sum(w * s for w, s in zip(weights, scores))
        
        return min(confidence, 1.0)
    
    def _compute_read_support_score(self, read_support: int) -> float:
        """Compute read support score."""
        if read_support < self.params.min_reads:
            return 0.0
        
        # Logarithmic scaling for read support
        import math
        score = min(math.log(read_support + 1) / math.log(10), 1.0)
        
        return score
    
    def _compute_coverage_score(self, coverage: float) -> float:
        """Compute coverage score."""
        if coverage < self.params.min_coverage:
            return 0.0
        
        return min(coverage, 1.0)
    
    def _compute_quality_score(self, quality_scores: List[float]) -> float:
        """Compute quality score."""
        if not quality_scores:
            return 0.0
        
        avg_quality = np.mean(quality_scores)
        
        if avg_quality < self.params.min_quality:
            return 0.0
        
        return min(avg_quality, 1.0)
    
    def _compute_consensus_score(self, splice_consensus: float) -> float:
        """Compute splice site consensus score."""
        return min(splice_consensus, 1.0)
    
    def _compute_splice_consensus(self, isoform: Isoform) -> float:
        """Compute splice site consensus score."""
        if not isoform.splice_sites:
            return 1.0  # Mono-exon isoforms
        
        consensus_scores = []
        
        for site in isoform.splice_sites:
            # Count supporting reads for this site
            supporting_reads = 0
            
            for read in isoform.read_signatures:
                for read_site in read.splice_sites:
                    if (abs(read_site.position - site.position) <= 10 and
                        read_site.donor_acceptor == site.donor_acceptor):
                        supporting_reads += 1
                        break
            
            # Compute consensus score for this site
            consensus_score = supporting_reads / len(isoform.read_signatures)
            consensus_scores.append(consensus_score)
        
        return np.mean(consensus_scores) if consensus_scores else 0.0


class AdvancedScorer:
    """Advanced confidence scoring with additional metrics."""
    
    def __init__(self, params: ScoringParams):
        self.params = params
        self.base_scorer = ConfidenceScorer(params)
    
    def score_isoform(self, isoform: Isoform) -> Dict[str, float]:
        """
        Compute comprehensive confidence scores for an isoform.
        
        Args:
            isoform: Isoform to score
        
        Returns:
            Dictionary of scores and metrics
        """
        # Base confidence score
        base_score = self.base_scorer.score_isoform(isoform)
        
        # Additional metrics
        metrics = {
            "base_confidence": base_score,
            "read_support": len(isoform.read_signatures),
            "coverage": isoform.coverage,
            "avg_quality": np.mean([r.quality_score for r in isoform.read_signatures]) if isoform.read_signatures else 0.0,
            "splice_consensus": self._compute_splice_consensus(isoform),
            "position_consensus": self._compute_position_consensus(isoform),
            "length_consensus": self._compute_length_consensus(isoform),
            "strand_consensus": self._compute_strand_consensus(isoform)
        }
        
        # Compute composite scores
        metrics["composite_score"] = self._compute_composite_score(metrics)
        metrics["reliability_score"] = self._compute_reliability_score(metrics)
        
        return metrics
    
    def _compute_splice_consensus(self, isoform: Isoform) -> float:
        """Compute splice site consensus."""
        if not isoform.splice_sites:
            return 1.0
        
        consensus_scores = []
        
        for site in isoform.splice_sites:
            supporting_reads = 0
            
            for read in isoform.read_signatures:
                for read_site in read.splice_sites:
                    if (abs(read_site.position - site.position) <= 10 and
                        read_site.donor_acceptor == site.donor_acceptor):
                        supporting_reads += 1
                        break
            
            consensus_score = supporting_reads / len(isoform.read_signatures)
            consensus_scores.append(consensus_score)
        
        return np.mean(consensus_scores) if consensus_scores else 0.0
    
    def _compute_position_consensus(self, isoform: Isoform) -> float:
        """Compute position consensus."""
        if not isoform.read_signatures:
            return 0.0
        
        starts = [r.start for r in isoform.read_signatures]
        ends = [r.end for r in isoform.read_signatures]
        
        # Compute coefficient of variation
        start_cv = np.std(starts) / np.mean(starts) if np.mean(starts) > 0 else 1.0
        end_cv = np.std(ends) / np.mean(ends) if np.mean(ends) > 0 else 1.0
        
        # Convert to consensus score (lower CV = higher consensus)
        position_consensus = 1.0 / (1.0 + start_cv + end_cv)
        
        return min(position_consensus, 1.0)
    
    def _compute_length_consensus(self, isoform: Isoform) -> float:
        """Compute length consensus."""
        if not isoform.read_signatures:
            return 0.0
        
        lengths = [r.end - r.start for r in isoform.read_signatures]
        
        # Compute coefficient of variation
        length_cv = np.std(lengths) / np.mean(lengths) if np.mean(lengths) > 0 else 1.0
        
        # Convert to consensus score
        length_consensus = 1.0 / (1.0 + length_cv)
        
        return min(length_consensus, 1.0)
    
    def _compute_strand_consensus(self, isoform: Isoform) -> float:
        """Compute strand consensus."""
        if not isoform.read_signatures:
            return 0.0
        
        strands = [r.strand for r in isoform.read_signatures]
        
        # Count strand agreement
        positive_strand = strands.count('+')
        negative_strand = strands.count('-')
        
        total_reads = len(strands)
        max_strand_count = max(positive_strand, negative_strand)
        
        strand_consensus = max_strand_count / total_reads
        
        return strand_consensus
    
    def _compute_composite_score(self, metrics: Dict[str, float]) -> float:
        """Compute composite score from all metrics."""
        weights = {
            "base_confidence": 0.3,
            "splice_consensus": 0.2,
            "position_consensus": 0.15,
            "length_consensus": 0.15,
            "strand_consensus": 0.1,
            "avg_quality": 0.1
        }
        
        composite_score = sum(weights.get(key, 0) * value for key, value in metrics.items())
        
        return min(composite_score, 1.0)
    
    def _compute_reliability_score(self, metrics: Dict[str, float]) -> float:
        """Compute reliability score based on consistency."""
        # Factors that affect reliability
        reliability_factors = [
            metrics["splice_consensus"],
            metrics["position_consensus"],
            metrics["length_consensus"],
            metrics["strand_consensus"]
        ]
        
        # Reliability is the minimum of all factors
        reliability_score = min(reliability_factors)
        
        return reliability_score


def score_isoforms(
    isoforms: List[Isoform],
    params: Optional[ScoringParams] = None,
    advanced: bool = False
) -> List[Dict[str, float]]:
    """
    Score a list of isoforms.
    
    Args:
        isoforms: List of isoforms to score
        params: Scoring parameters (uses defaults if None)
        advanced: Whether to use advanced scoring
    
    Returns:
        List of score dictionaries
    """
    if params is None:
        params = ScoringParams()
    
    if advanced:
        scorer = AdvancedScorer(params)
        return [scorer.score_isoform(isoform) for isoform in isoforms]
    else:
        scorer = ConfidenceScorer(params)
        return [{"confidence": scorer.score_isoform(isoform)} for isoform in isoforms]


def filter_by_confidence(
    isoforms: List[Isoform],
    min_confidence: float = 0.5,
    params: Optional[ScoringParams] = None
) -> List[Isoform]:
    """
    Filter isoforms by confidence score.
    
    Args:
        isoforms: List of isoforms to filter
        min_confidence: Minimum confidence threshold
        params: Scoring parameters (uses defaults if None)
    
    Returns:
        List of filtered isoforms
    """
    if params is None:
        params = ScoringParams()
    
    scorer = ConfidenceScorer(params)
    
    filtered_isoforms = []
    for isoform in isoforms:
        confidence = scorer.score_isoform(isoform)
        if confidence >= min_confidence:
            isoform.confidence_score = confidence
            filtered_isoforms.append(isoform)
    
    return filtered_isoforms
