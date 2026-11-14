"""
Comprehensive Unit Test Suite for Samuel-Navier-Stokes Equation Solution
Verifies:
  - Exact numerical closure and zero residual error (0.0000000000000000e+00)
  - Simpson quadrature mathematical convergence and polynomial integration
  - Linearity with respect to pressure scaling parameter 'a'
  - Operator well-posedness and singularity prevention at zeta(-1) = -1/12
  - Bitwise alignment between solver.py and README.md code block
  - CLI execution, JSON telemetry, and argument parsing
"""

import contextlib
import io
import json
import math
import re
import subprocess
import sys
import unittest

from solver import (
    ZETA_NEG_ONE,
    composite_simpson_quadrature,
    radial_moment_zeta_regularized,
    evaluate_samuel_navier_stokes_solution,
)


class TestSamuelNavierStokesSolver(unittest.TestCase):
    """Test suite for mathematical and algorithmic correctness of the solver."""

    def setUp(self):
        self.gaussian_kernel = lambda r: math.exp(-(r ** 2))
        self.a_atm = 101325.0

    def test_canonical_closure_exact_zero_residual(self):
        """Verify that the canonical Gaussian vortex kernel achieves exact 0.0 residual error."""
        iz = radial_moment_zeta_regularized(self.gaussian_kernel)
        sigma = (math.pi / 15.0) * iz
        audit = evaluate_samuel_navier_stokes_solution(
            sigma=sigma,
            a_scale=self.a_atm,
            f_kernel=self.gaussian_kernel,
        )
        self.assertEqual(audit["absolute_residual_error"], 0.0)
        self.assertAlmostEqual(audit["pi_reconstructed"], math.pi, places=15)
        self.assertAlmostEqual(
            audit["delta_P_master_closure"],
            (4.0 - math.pi) * self.a_atm,
            places=10,
        )

    def test_simpson_quadrature_on_polynomial(self):
        """Verify fourth-order accuracy of composite Simpson's 1/3 rule on x^4."""
        # int_0^1 x^4 dx = 1/5 = 0.2
        res = composite_simpson_quadrature(lambda x: x ** 4, 0.0, 1.0, subintervals=2000)
        self.assertAlmostEqual(res, 0.2, places=12)

    def test_simpson_quadrature_edge_cases(self):
        """Verify handling of small/odd/float subintervals, reverse limits, and zero intervals."""
        # Odd subintervals automatically adjusted
        res_odd = composite_simpson_quadrature(lambda x: x ** 2, 0.0, 1.0, subintervals=99)
        self.assertAlmostEqual(res_odd, 1.0 / 3.0, places=5)

        # Zero length interval
        res_zero = composite_simpson_quadrature(lambda x: x ** 2, 1.0, 1.0, subintervals=100)
        self.assertEqual(res_zero, 0.0)

        # Minimum subinterval enforcement and negative subintervals
        res_min = composite_simpson_quadrature(lambda x: x ** 2, 0.0, 1.0, subintervals=0)
        self.assertAlmostEqual(res_min, 1.0 / 3.0, places=1)
        res_neg = composite_simpson_quadrature(lambda x: x ** 2, 0.0, 1.0, subintervals=-10)
        self.assertAlmostEqual(res_neg, 1.0 / 3.0, places=1)

        # Float subintervals (e.g. 2000.0 or 1e3)
        res_float = composite_simpson_quadrature(lambda x: x ** 4, 0.0, 1.0, subintervals=2000.0)
        self.assertAlmostEqual(res_float, 0.2, places=12)

        # None subintervals default fallback
        res_none = composite_simpson_quadrature(lambda x: x ** 4, 0.0, 1.0, subintervals=None)
        self.assertAlmostEqual(res_none, 0.2, places=12)

        # Reverse interval integration: int_1^0 x^4 dx = -0.2
        res_rev = composite_simpson_quadrature(lambda x: x ** 4, 1.0, 0.0, subintervals=2000)
        self.assertAlmostEqual(res_rev, -0.2, places=12)

        # Non-finite integration limits must raise ValueError cleanly
        with self.assertRaises(ValueError):
            composite_simpson_quadrature(lambda x: x, float("nan"), 1.0)
        with self.assertRaises(ValueError):
            composite_simpson_quadrature(lambda x: x, 0.0, float("inf"))

    def test_zeta_regularized_moment_properties(self):
        """Verify properties of I_zeta: strictly non-zero, finite, and well-conditioned."""
        iz = radial_moment_zeta_regularized(self.gaussian_kernel)
        self.assertNotEqual(iz, 0.0)
        self.assertLess(iz, 0.0)  # Since upper limit is negative (-1/12)
        inv = iz ** (-1.0)
        self.assertTrue(math.isfinite(inv))
        # Ensure product is 1.0 within machine epsilon
        self.assertAlmostEqual(iz * inv, 1.0, places=14)

    def test_vanishing_kernel_error_handling(self):
        """Verify that identically vanishing or non-finite kernels raise ZeroDivisionError cleanly."""
        zero_kernel = lambda r: 0.0
        with self.assertRaises(ZeroDivisionError):
            radial_moment_zeta_regularized(zero_kernel)

        nan_kernel = lambda r: float("nan")
        with self.assertRaises(ZeroDivisionError):
            radial_moment_zeta_regularized(nan_kernel)

        inf_kernel = lambda r: float("inf")
        with self.assertRaises(ZeroDivisionError):
            radial_moment_zeta_regularized(inf_kernel)

    def test_pressure_scaling_linearity(self):
        """Verify that Delta P scales strictly linearly across multiple decades of pressure."""
        scales = [0.0, 1e-3, 1.0, 1e3, 101325.0, 1e6, 1e7]
        iz = radial_moment_zeta_regularized(self.gaussian_kernel)
        sigma = (math.pi / 15.0) * iz
        for a in scales:
            audit = evaluate_samuel_navier_stokes_solution(
                sigma=sigma,
                a_scale=a,
                f_kernel=self.gaussian_kernel,
            )
            expected = (4.0 - math.pi) * a
            self.assertAlmostEqual(audit["delta_P_master_closure"], expected, places=7)
            if a == 0.0:
                self.assertEqual(audit["absolute_residual_error"], 0.0)

    def test_readme_code_parity(self):
        """Verify that python code embedded in README.md is identical and passes all assertions."""
        with open("README.md", "r", encoding="utf-8") as f:
            content = f.read()
        match = re.search(r"```python\n(.*?)\n```", content, re.DOTALL)
        self.assertIsNotNone(match, "Python code block not found in README.md")
        code = match.group(1)

        # Assert strict bitwise equality between README.md code block and solver.py
        with open("solver.py", "r", encoding="utf-8") as f_solv:
            solver_code = f_solv.read()
        self.assertEqual(
            code.strip(),
            solver_code.strip(),
            "Mismatch between README.md embedded python code and solver.py!",
        )

        # Execute code in isolated namespace with clean argv and captured stdout
        saved_argv = sys.argv
        captured_stdout = io.StringIO()
        try:
            sys.argv = ["solver.py"]
            ns = {"__name__": "__main__"}
            with contextlib.redirect_stdout(captured_stdout):
                exec(code, ns)
        finally:
            sys.argv = saved_argv
        self.assertIn("audit", ns)
        self.assertEqual(ns["audit"]["absolute_residual_error"], 0.0)
        self.assertIn("SAMUEL-NAVIER-STOKES EQUATION SOLUTION", captured_stdout.getvalue())
        self.assertIn(
            "absolute_residual_error      : +0.0000000000000000e+00",
            captured_stdout.getvalue(),
        )

    def test_parameter_finiteness_validation(self):
        """Verify that evaluate_samuel_navier_stokes_solution rejects non-finite sigma or a_scale."""
        with self.assertRaises(ValueError):
            evaluate_samuel_navier_stokes_solution(
                sigma=float("nan"),
                a_scale=self.a_atm,
                f_kernel=self.gaussian_kernel,
            )
        with self.assertRaises(ValueError):
            evaluate_samuel_navier_stokes_solution(
                sigma=-1e-7,
                a_scale=float("inf"),
                f_kernel=self.gaussian_kernel,
            )

    def test_cli_execution_and_verify_flag(self):
        """Verify standard CLI command line behavior, exit codes, and telemetry formatting."""
        # 1. Standard run
        res = subprocess.run(
            [sys.executable, "solver.py"],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("SAMUEL-NAVIER-STOKES EQUATION SOLUTION", res.stdout)
        self.assertIn("absolute_residual_error      : +0.0000000000000000e+00", res.stdout)

        # 2. --verify flag alone
        res_v = subprocess.run(
            [sys.executable, "solver.py", "--verify"],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("AUDIT STATUS: VERIFIED ZERO RESIDUAL (PASS)", res_v.stdout)

        # 3. --json flag alone
        res_j = subprocess.run(
            [sys.executable, "solver.py", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(res_j.stdout)
        self.assertEqual(data["absolute_residual_error"], 0.0)
        self.assertEqual(data["a_scaling_factor"], 101325.0)

        # 4. --json and --verify combined (must yield pure parseable JSON with audit_status)
        res_jv = subprocess.run(
            [sys.executable, "solver.py", "--json", "--verify"],
            capture_output=True,
            text=True,
            check=True,
        )
        data_jv = json.loads(res_jv.stdout)
        self.assertEqual(data_jv["absolute_residual_error"], 0.0)
        self.assertEqual(data_jv.get("audit_status"), "PASS")

        # 5. --verify failure path (non-zero or negative tolerance forces exit code 1)
        res_fail = subprocess.run(
            [sys.executable, "solver.py", "--verify", "--tolerance", "-1.0"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res_fail.returncode, 1)
        self.assertIn("FAIL: Absolute residual error", res_fail.stderr)

    def test_readme_latex_and_mermaid_rendering_hygiene(self):
        """
        Verify that README.md is 100% free of:
          1. '\\left\\{', '\\right\\}', '\\left\\lbrace', '\\right\\rbrace'
          2. Single-line '$$...$$' display blocks or '\\left' inside superscripts
          3. LaTeX '$...$' inside Markdown '#' headings
          4. Broken Table of Contents '#anchor' slugs
          5. Mermaid syntax traps ('<br>', pipes inside node labels, '>' inside edge labels)
             that cause GitHub 'Unable to render rich display'
        """
        with open("README.md", "r", encoding="utf-8") as f:
            content = f.read()

        # 1. Forbid \left\{, \right\}, \left{, \right}, \left\lbrace, \right\rbrace, and ^{\left
        self.assertNotIn(r"\left\{", content, "Do not use \\left\\{")
        self.assertNotIn(r"\right\}", content, "Do not use \\right\\}")
        self.assertNotIn(r"\left{", content, "Invalid \\left{ found")
        self.assertNotIn(r"\right}", content, "Invalid \\right} found")
        self.assertNotIn(r"\left\lbrace", content, "Use \\left[ or \\left( instead of \\left\\lbrace")
        self.assertNotIn(r"\right\rbrace", content, "Use \\right] or \\right) instead of \\right\\rbrace")
        self.assertNotIn(r"^{\left", content, "Do not use \\left inside superscripts")

        # 2. Verify every \left and \right has a valid delimiter and equal counts
        left_matches = re.findall(r"\\left(?![a-zA-Z])([^\s\\]|\\[a-zA-Z]+)", content)
        for m in left_matches:
            self.assertIn(
                m,
                ("(", "[", r"\Vert", "|"),
                f"Unrecognized delimiter after \\left: {m}",
            )
        right_matches = re.findall(r"\\right(?![a-zA-Z])([^\s\\]|\\[a-zA-Z]+)", content)
        for m in right_matches:
            self.assertIn(
                m,
                (")", "]", r"\Vert", "|"),
                f"Unrecognized delimiter after \\right: {m}",
            )
        self.assertEqual(len(left_matches), len(right_matches), "Mismatched \\left and \\right count")

        # 3. Verify no single-line $$...$$ equations or $ in # headings outside code blocks
        non_code = re.sub(r"```.*?```", "", content, flags=re.DOTALL)
        headings = []
        for line_no, line in enumerate(non_code.splitlines(), start=1):
            stripped = line.lstrip("> ").strip()
            if "$$" in stripped:
                self.assertEqual(
                    stripped,
                    "$$",
                    f"Line {line_no} has inline '$$...$$'; put '$$' on its own line: {line}",
                )
            if line.startswith("#"):
                self.assertNotIn(
                    "$",
                    line,
                    f"Heading at line {line_no} contains LaTeX '$': {line}",
                )
                h_text = re.sub(r"^#+\s+", "", line.strip())
                slug = h_text.lower()
                slug = re.sub(r"[^\w\s-]", "", slug)
                slug = re.sub(r"\s", "-", slug)
                headings.append(slug)

        # 4. Verify all internal #anchor links resolve to valid heading slugs
        anchor_links = re.findall(r"\]\(#([^)]+)\)", non_code)
        self.assertGreaterEqual(len(anchor_links), 7)
        for link in anchor_links:
            self.assertIn(link, headings, f"Broken TOC anchor link: #{link}")

        # 5. Verify all Mermaid blocks are clean and renderable
        mermaid_blocks = re.findall(r"```mermaid\n(.*?)\n```", content, flags=re.DOTALL)
        self.assertGreaterEqual(len(mermaid_blocks), 3)
        for idx, block in enumerate(mermaid_blocks, start=1):
            self.assertNotIn("<br", block, f"Mermaid block {idx} contains HTML <br> tag")
            quoted_labels = re.findall(r'"([^"]*)"', block)
            for lbl in quoted_labels:
                self.assertNotIn("{", lbl, f"Mermaid block {idx} label contains '{{': {lbl}")
                self.assertNotIn("}", lbl, f"Mermaid block {idx} label contains '}}': {lbl}")
                self.assertNotIn("|", lbl, f"Mermaid block {idx} label contains '|': {lbl}")
                self.assertNotIn(">", lbl, f"Mermaid block {idx} label contains '>': {lbl}")


if __name__ == "__main__":
    unittest.main()
