"""v2: protein sequences, CGC clusters, ARG + GO annotations

Adds the BioForge v2 schema (all additive):
  - genes.protein_seq, genes.cluster_id
  - cazyme_clusters  (dbCAN CGC + substrate)
  - arg_annotations  (hAMRonization AMR genes)
  - go_annotations   (normalised GO terms)

Revision ID: b2c3d4e5f6a7
Revises: 067f9f8305f6
Create Date: 2026-07-07
"""
from alembic import op
import sqlalchemy as sa


revision = 'b2c3d4e5f6a7'
down_revision = '067f9f8305f6'
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        'cazyme_clusters',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('sample_id', sa.Integer(), nullable=False),
        sa.Column('cluster_key', sa.String(length=64), nullable=False),
        sa.Column('contig', sa.String(length=255), nullable=False),
        sa.Column('start', sa.Integer(), nullable=False),
        sa.Column('end', sa.Integer(), nullable=False),
        sa.Column('composition', sa.String(length=128), nullable=True),
        sa.Column('n_genes', sa.Integer(), nullable=False),
        sa.Column('substrate', sa.Text(), nullable=True),
        sa.Column('substrate_score', sa.Float(), nullable=True),
        sa.Column('raw', sa.JSON(), nullable=False),
        sa.ForeignKeyConstraint(['sample_id'], ['samples.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
    )
    with op.batch_alter_table('cazyme_clusters', schema=None) as batch_op:
        batch_op.create_index('ix_cgc_key', ['cluster_key'], unique=False)
        batch_op.create_index('ix_cgc_locus', ['sample_id', 'contig', 'start', 'end'],
                              unique=False)

    op.create_table(
        'arg_annotations',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('gene_id', sa.Integer(), nullable=True),
        sa.Column('sample_id', sa.Integer(), nullable=False),
        sa.Column('gene_symbol', sa.String(length=128), nullable=True),
        sa.Column('gene_name', sa.Text(), nullable=True),
        sa.Column('drug_class', sa.String(length=255), nullable=True),
        sa.Column('resistance_mechanism', sa.String(length=255), nullable=True),
        sa.Column('identity', sa.Float(), nullable=True),
        sa.Column('coverage', sa.Float(), nullable=True),
        sa.Column('tool', sa.String(length=64), nullable=True),
        sa.Column('reference_db', sa.String(length=64), nullable=True),
        sa.Column('accession', sa.String(length=64), nullable=True),
        sa.Column('raw', sa.JSON(), nullable=False),
        sa.ForeignKeyConstraint(['gene_id'], ['genes.id'], ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['sample_id'], ['samples.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
    )
    with op.batch_alter_table('arg_annotations', schema=None) as batch_op:
        batch_op.create_index('ix_arg_drug_class', ['drug_class'], unique=False)
        batch_op.create_index('ix_arg_symbol', ['gene_symbol'], unique=False)

    op.create_table(
        'go_annotations',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('gene_id', sa.Integer(), nullable=False),
        sa.Column('go_id', sa.String(length=20), nullable=False),
        sa.Column('source', sa.String(length=64), nullable=True),
        sa.ForeignKeyConstraint(['gene_id'], ['genes.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
    )
    with op.batch_alter_table('go_annotations', schema=None) as batch_op:
        batch_op.create_index('ix_go_id', ['go_id'], unique=False)

    with op.batch_alter_table('genes', schema=None) as batch_op:
        batch_op.add_column(sa.Column('protein_seq', sa.Text(), nullable=True))
        batch_op.add_column(sa.Column('cluster_id', sa.Integer(), nullable=True))
        batch_op.create_foreign_key(
            'fk_genes_cluster', 'cazyme_clusters', ['cluster_id'], ['id'],
            ondelete='SET NULL',
        )


def downgrade() -> None:
    with op.batch_alter_table('genes', schema=None) as batch_op:
        batch_op.drop_constraint('fk_genes_cluster', type_='foreignkey')
        batch_op.drop_column('cluster_id')
        batch_op.drop_column('protein_seq')

    with op.batch_alter_table('go_annotations', schema=None) as batch_op:
        batch_op.drop_index('ix_go_id')
    op.drop_table('go_annotations')

    with op.batch_alter_table('arg_annotations', schema=None) as batch_op:
        batch_op.drop_index('ix_arg_symbol')
        batch_op.drop_index('ix_arg_drug_class')
    op.drop_table('arg_annotations')

    with op.batch_alter_table('cazyme_clusters', schema=None) as batch_op:
        batch_op.drop_index('ix_cgc_locus')
        batch_op.drop_index('ix_cgc_key')
    op.drop_table('cazyme_clusters')
