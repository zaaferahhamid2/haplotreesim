# Contributing to HaploTreeSim

Thank you for your interest in contributing to HaploTreeSim! This document provides guidelines for contributing to the project.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/haplotreesim.git
   cd haplotreesim
   ```
3. **Create a branch** for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Setup

1. Install in development mode:
   ```bash
   pip install -e ".[dev]"
   ```

2. Run tests to ensure everything works:
   ```bash
   python tests/test_week5.py
   ```

## Making Changes

### Code Style

- Follow [PEP 8](https://pep8.org/) style guidelines
- Use type hints for function signatures
- Add docstrings to all public functions and classes
- Keep lines under 100 characters when possible
- Use descriptive variable names

Example:
```python
def calculate_copy_number(
    clone: Clone,
    bin_index: int
) -> Tuple[int, int]:
    """
    Calculate haplotype-specific copy numbers for a clone at a bin.
    
    Args:
        clone: Clone object containing CN profiles
        bin_index: Index of the genomic bin
        
    Returns:
        Tuple of (cn_A, cn_B) copy numbers
    """
    return clone.cn_profile_A[bin_index], clone.cn_profile_B[bin_index]
```

### Testing

- Add tests for all new features
- Ensure all existing tests pass
- Aim for >80% code coverage
- Test edge cases and error conditions

Run tests:
```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=haplotreesim --cov-report=html
```

### Documentation

- Update README.md if adding user-facing features
- Add docstrings to new functions/classes
- Update relevant documentation files in `docs/`
- Include usage examples for new features

## Submitting Changes

1. **Commit your changes** with clear, descriptive messages:
   ```bash
   git commit -m "Add feature: haplotype-aware event generation"
   ```

2. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Open a Pull Request** on GitHub:
   - Provide a clear title and description
   - Reference any related issues
   - Describe what testing you've done

### Pull Request Checklist

- [ ] Code follows project style guidelines
- [ ] All tests pass
- [ ] New tests added for new features
- [ ] Documentation updated
- [ ] Commit messages are clear and descriptive
- [ ] No merge conflicts

## Reporting Issues

When reporting bugs or requesting features:

1. **Search existing issues** to avoid duplicates
2. **Use issue templates** when available
3. **Provide details**:
   - For bugs: steps to reproduce, expected vs actual behavior, system info
   - For features: use case, proposed solution, alternatives considered

## Development Workflow

### Week-by-Week Implementation

This project follows a 28-week implementation plan. Check the current week's focus in the README or project board.

### Branch Naming

- `feature/` - New features (e.g., `feature/week6-cna-events`)
- `bugfix/` - Bug fixes (e.g., `bugfix/negative-binomial-overflow`)
- `docs/` - Documentation updates
- `test/` - Test additions/improvements

### Commit Messages

Use clear, imperative commit messages:

✅ Good:
- `Add WGD event generation to simulator`
- `Fix allele count overflow in Beta-Binomial model`
- `Update README with installation instructions`

❌ Bad:
- `Update code`
- `Fixed bug`
- `Changes`

## Code Review Process

1. Maintainers will review your PR
2. Address feedback by pushing new commits
3. Once approved, your PR will be merged
4. Your contribution will be credited in the changelog

## Questions?

- Open an issue for general questions
- Tag issues with `question` label
- Check existing issues and documentation first

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to HaploTreeSim! 🧬
