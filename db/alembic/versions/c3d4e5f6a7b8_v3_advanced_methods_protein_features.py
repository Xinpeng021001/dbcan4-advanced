"""v3: advanced CAZyme methods + protein features

Adds the dbCAN4-advanced integration schema (all additive; baseline rows keep
working unchanged):
  - cazyme_annotations.confidence      (Float, nullable)   calibrated 0–1 score
  - cazyme_annotations.method_family   (String, 'baseline'|'advanced')
  - cazyme_annotations.method_kind     (String, nullable)  hmm/…/sequence-plm/…
  - cazyme_annotations.release_id      (FK releases, nullable)
  - protein_features                   (new table: SignalP6 / DeepTMHMM / structure)

Existing baseline rows are backfilled to method_family='baseline' via the column
server_default; the ingester also sets release_id/method_kind going forward.

Revision ID: c3d4e5f6a7b8
Revises: b2c3d4e5f6a7
Create Date: 2026-07-08
"""
from alembic import op
import sqlalchemy as sa


revision = 'c3d4e5f6a7b8'
down_revision = 'b2c3d4e5f6a7'
branch_labels = None
depends_on = None


def upgrade() -> None:
    # --- cazyme_annotations: advanced-method columns (additive) ---
    with op.batch_alter_table('cazyme_annotations', schema=None) as batch_op:
        batch_op.add_column(sa.Column('confidence', sa.Float(), nullable=True))
        batch_op.add_column(sa.Column(
            'method_family', sa.String(length=16),
            nullable=False, server_default='baseline'))
        batch_op.add_column(sa.Column('method_kind', sa.String(length=24), nullable=True))
        batch_op.add_column(sa.Column('release_id', sa.Integer(), nullable=True))
        batch_op.create_index('ix_cazyme_method_family', ['method_family'], unique=False)
        batch_op.create_index('ix_cazyme_release', ['release_id'], unique=False)
        batch_op.create_foreign_key(
            'fk_cazyme_release', 'releases', ['release_id'], ['id'],
            ondelete='CASCADE')

    # --- protein_features: SignalP6 / DeepTMHMM / structure / … ---
    op.create_table(
        'protein_features',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('gene_id', sa.Integer(), nullable=False),
        sa.Column('feature_type', sa.String(length=32), nullable=False),
        sa.Column('tool', sa.String(length=48), nullable=True),
        sa.Column('label', sa.String(length=128), nullable=True),
        sa.Column('score', sa.Float(), nullable=True),
        sa.Column('start', sa.Integer(), nullable=True),
        sa.Column('end', sa.Integer(), nullable=True),
        sa.Column('structure_path', sa.Text(), nullable=True),
        sa.Column('attributes', sa.JSON(), nullable=False),
        sa.Column('release_id', sa.Integer(), nullable=True),
        sa.ForeignKeyConstraint(['gene_id'], ['genes.id'], ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['release_id'], ['releases.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
    )
    with op.batch_alter_table('protein_features', schema=None) as batch_op:
        batch_op.create_index('ix_pf_gene', ['gene_id'], unique=False)
        batch_op.create_index('ix_pf_type', ['feature_type'], unique=False)
        batch_op.create_index('ix_pf_release', ['release_id'], unique=False)


def downgrade() -> None:
    with op.batch_alter_table('protein_features', schema=None) as batch_op:
        batch_op.drop_index('ix_pf_release')
        batch_op.drop_index('ix_pf_type')
        batch_op.drop_index('ix_pf_gene')
    op.drop_table('protein_features')

    with op.batch_alter_table('cazyme_annotations', schema=None) as batch_op:
        batch_op.drop_constraint('fk_cazyme_release', type_='foreignkey')
        batch_op.drop_index('ix_cazyme_release')
        batch_op.drop_index('ix_cazyme_method_family')
        batch_op.drop_column('release_id')
        batch_op.drop_column('method_kind')
        batch_op.drop_column('method_family')
        batch_op.drop_column('confidence')
